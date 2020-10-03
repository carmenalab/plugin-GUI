/*
    ------------------------------------------------------------------

    This file is part of the Open Ephys GUI
    Copyright (C) 2016 Open Ephys

    ------------------------------------------------------------------

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <stdio.h>

#include "CAR.h"
#include "CAREditor.h"
#include <algorithm>



CAR::CAR()
    : GenericProcessor ("Common Avg Ref") //, threshold(200.0), state(true)
{
    setProcessorType (PROCESSOR_TYPE_FILTER);

    m_avgBuffer = AudioSampleBuffer (1, 10000); // 1-dimensional buffer to hold the avg
    m_mode = Mode::MEAN;
}


CAR::~CAR()
{
}


AudioProcessorEditor* CAR::createEditor()
{
    editor = new CAREditor (this, true);
    return editor;
}


float CAR::getGainLevel()
{
    m_gainLevel.updateTarget();
    return m_gainLevel.getNextValue();
}


void CAR::setGainLevel (float newGain)
{
    m_gainLevel.setValue (newGain);
}


void CAR::process (AudioSampleBuffer& buffer)
{
    const int numSamples            = buffer.getNumSamples();
    const int numReferenceChannels  = m_referenceChannels.size();
    const int numAffectedChannels   = m_affectedChannels.size();

    // There are no sense to do any processing if either number of reference or affected channels is zero.
    if (! numReferenceChannels || ! numAffectedChannels)
    {
        return;
    }

    switch (m_mode) {
        case Mode::MEAN:
        default:
            computeSubtractorForMean(buffer);
            break;
        case Mode::MEDIAN:
            computeSubtractorForMedian(buffer);
            break;
    }

    m_gainLevel.updateTarget();
    const float gain =  -1.0f * m_gainLevel.getNextValue() / 100.f;

    for (int i = 0; i < numAffectedChannels; ++i)
    {
        buffer.addFrom (m_affectedChannels[i],  // destChannel
                        0,                      // destStartSample
                        m_avgBuffer,            // source
                        0,                      // sourceChannel
                        0,                      // sourceStartSample
                        numSamples,             // numSamples
                        gain);                  // gain to apply
    }
}


void CAR::setReferenceChannels (const Array<int>& newReferenceChannels)
{
    const ScopedLock myScopedLock (objectLock);

    m_referenceChannels = Array<int> (newReferenceChannels);
}


void CAR::setAffectedChannels (const Array<int>& newAffectedChannels)
{
    const ScopedLock myScopedLock (objectLock);

    m_affectedChannels = Array<int> (newAffectedChannels);
}


void CAR::setReferenceChannelState (int channel, bool newState)
{
    if (! newState)
        m_referenceChannels.removeFirstMatchingValue (channel);
    else
        m_referenceChannels.addIfNotAlreadyThere (channel);
}


void CAR::setAffectedChannelState (int channel, bool newState)
{
    if (! newState)
        m_affectedChannels.removeFirstMatchingValue (channel);
    else
        m_affectedChannels.addIfNotAlreadyThere (channel);
}

void CAR::saveCustomChannelParametersToXml(XmlElement* channelElement,
    int channelNumber, InfoObjectCommon::InfoObjectType channelType)
{
    if (channelType == InfoObjectCommon::DATA_CHANNEL)
    {
        XmlElement* groupState = channelElement->createNewChildElement("GROUPSTATE");

        const Array<int>& referenceChannels = getReferenceChannels();
        bool isReferenceChannel = referenceChannels.contains(channelNumber);
        groupState->setAttribute("reference", isReferenceChannel);

        const Array<int>& affectedChannels = getAffectedChannels();
        bool isAffectedChannel = affectedChannels.contains(channelNumber);
        groupState->setAttribute("affected", isAffectedChannel);
    }
}

void CAR::loadCustomChannelParametersFromXml(XmlElement* channelElement,
    InfoObjectCommon::InfoObjectType channelType)
{
    if (channelType == InfoObjectCommon::DATA_CHANNEL)
    {
        int channelNumber = channelElement->getIntAttribute("number");

        forEachXmlChildElementWithTagName(*channelElement, groupState, "GROUPSTATE")
        {
            if (groupState->hasAttribute("reference"))
            {
                bool isReferenceChannel = groupState->getBoolAttribute("reference");
                setReferenceChannelState(channelNumber, isReferenceChannel);
            }

            if (groupState->hasAttribute("affected"))
            {
                bool isAffectedChannel = groupState->getBoolAttribute("affected");
                setAffectedChannelState(channelNumber, isAffectedChannel);
            }
        }
    }
}

void CAR::computeSubtractorForMean(const AudioSampleBuffer &buffer) {
    const int numSamples = buffer.getNumSamples();
    const int numReferenceChannels  = m_referenceChannels.size();

    m_avgBuffer.clear();

    for (int i = 0; i < numReferenceChannels; ++i)
    {
        m_avgBuffer.addFrom (0,                         // destChannel
                             0,                         // destStartSample
                             buffer,                    // source
                             m_referenceChannels[i],    // sourceChannel
                             0,                         // sourceStartSample
                             numSamples,                // numSamples
                             1.0f);                     // gain to apply
    }

    m_avgBuffer.applyGain (1.0f / float (numReferenceChannels));
}

void CAR::computeSubtractorForMedian(const AudioSampleBuffer &buffer) {
    const int numSamples = buffer.getNumSamples();
    const int numReferenceChannels  = m_referenceChannels.size();

    m_avgBuffer.clear();

    m_medianBuffer.resize(numReferenceChannels);

    for (int sampleIdx = 0; sampleIdx < numSamples; sampleIdx++) {
        // For each sample index, copy in each channel's data to a temporary buffer;
        // then sort that cross-channel data to get the median.
        for (int chidx = 0; chidx < numReferenceChannels; chidx++) {
            m_medianBuffer[chidx] = buffer.getSample(m_referenceChannels[chidx], sampleIdx);
        }
        std::sort(m_medianBuffer.begin(), m_medianBuffer.end());

        float median = m_medianBuffer[(int) numReferenceChannels / 2];
        if (numReferenceChannels % 2 == 0) {
            median = (median + m_medianBuffer[((int) numReferenceChannels / 2) - 1]) / 2.0f;
        }
        m_avgBuffer.setSample(0, sampleIdx, median);
    }
}
