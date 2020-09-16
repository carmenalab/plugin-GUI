/*
    ------------------------------------------------------------------

    This file is part of the Open Ephys GUI
    Copyright (C) 2014 Open Ephys

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

#include "FilterEditor.h"
#include "FilterNode.h"
#include <stdio.h>


FilterEditor::FilterEditor(GenericProcessor* parentNode, bool useDefaultParameterEditors=true)
    : GenericEditor(parentNode, useDefaultParameterEditors)

{
    desiredWidth = 150;

    lastLowCutString = " ";
    lastHighCutString = " ";
    lastOrderString = " ";

    highCutLabel = new Label("high cut label", "High cut:");
    highCutLabel->setBounds(10,65,80,20);
    highCutLabel->setFont(Font("Small Text", 12, Font::plain));
    highCutLabel->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(highCutLabel);

    lowCutLabel = new Label("low cut label", "Low cut:");
    lowCutLabel->setBounds(10,25,80,20);
    lowCutLabel->setFont(Font("Small Text", 12, Font::plain));
    lowCutLabel->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(lowCutLabel);

    lowCutValue = new Label("low cut value", lastLowCutString);
    lowCutValue->setBounds(15,42,60,18);
    lowCutValue->setFont(Font("Default", 15, Font::plain));
    lowCutValue->setColour(Label::textColourId, Colours::white);
    lowCutValue->setColour(Label::backgroundColourId, Colours::grey);
    lowCutValue->setEditable(true);
    lowCutValue->addListener(this);
    lowCutValue->setTooltip("Set the low cut for the selected channels");
    addAndMakeVisible(lowCutValue);

    highCutValue = new Label("high cut label", lastHighCutString);
    highCutValue->setBounds(15,82,60,18);
    highCutValue->setFont(Font("Default", 15, Font::plain));
    highCutValue->setColour(Label::textColourId, Colours::white);
    highCutValue->setColour(Label::backgroundColourId, Colours::grey);
    highCutValue->setEditable(true);
    highCutValue->addListener(this);
    highCutValue->setTooltip("Set the high cut for the selected channels");
    addAndMakeVisible(highCutValue);

    orderLabel = new Label("order label", "Order:");
    orderLabel->setBounds(85,25,80,20);
    orderLabel->setFont(Font("Small Text", 12, Font::plain));
    orderLabel->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(orderLabel);

    orderValue = new Label("filter order value", lastOrderString);
    orderValue->setBounds(95,42,30,18);
    orderValue->setFont(Font("Default", 15, Font::plain));
    orderValue->setColour(Label::textColourId, Colours::white);
    orderValue->setColour(Label::backgroundColourId, Colours::grey);
    orderValue->setEditable(true);
    orderValue->addListener(this);
    orderValue->setTooltip("Set the filter order for the selected channels");
    addAndMakeVisible(orderValue);


    applyFilterOnADC = new UtilityButton("+ADCs",Font("Default", 10, Font::plain));
    applyFilterOnADC->addListener(this);
    applyFilterOnADC->setBounds(90,70,40,18);
    applyFilterOnADC->setClickingTogglesState(true);
    applyFilterOnADC->setTooltip("When this button is off, ADC channels will not be filtered");
    addAndMakeVisible(applyFilterOnADC);

    applyFilterOnChan = new UtilityButton("+CH",Font("Default", 10, Font::plain));
    applyFilterOnChan->addListener(this);
    applyFilterOnChan->setBounds(95,95,30,18);
    applyFilterOnChan->setClickingTogglesState(true);
    applyFilterOnChan->setToggleState(true, dontSendNotification);
    applyFilterOnChan->setTooltip("When this button is off, selected channels will not be filtered");
    addAndMakeVisible(applyFilterOnChan);

}

FilterEditor::~FilterEditor()
{

}

void FilterEditor::setDefaults(double lowCut, double highCut, int order)
{
    lastHighCutString = String(roundFloatToInt(highCut));
    lastLowCutString = String(roundFloatToInt(lowCut));
    lastOrderString = String(order);

    resetToSavedText();
}

void FilterEditor::resetToSavedText()
{
    highCutValue->setText(lastHighCutString, dontSendNotification);
    lowCutValue->setText(lastLowCutString, dontSendNotification);
    orderValue->setText(lastOrderString, dontSendNotification);
}


void FilterEditor::labelTextChanged(Label* label)
{
    FilterNode* fn = (FilterNode*) getProcessor();

    Value val = label->getTextValue();
    double requestedValue = double(val.getValue());

    if (requestedValue < 0.01 || requestedValue > FilterNode::MAX_FREQUENCY)
    {
        CoreServices::sendStatusMessage("Value out of range.");

        if (label == highCutValue)
        {
            label->setText(lastHighCutString, dontSendNotification);
            lastHighCutString = label->getText();
        }
        else
        {
            label->setText(lastLowCutString, dontSendNotification);
            lastLowCutString = label->getText();
        }

        return;
    }

    Array<int> chans = getActiveChannels();

    // This needs to change, since there's not enough feedback about whether
    // or not individual channel settings were altered:

    for (int n = 0; n < chans.size(); n++)
    {

        if (label == highCutValue)
        {
            double minVal = fn->getLowCutValueForChannel(chans[n]);

            if (requestedValue > minVal)
            {
                fn->setCurrentChannel(chans[n]);
                fn->setParameter(FilterNode::PARAMETER_INDEX_HIGH_CUT, requestedValue);
            }

            lastHighCutString = label->getText();

        }
        else if (label == lowCutValue)
        {
            double maxVal = fn->getHighCutValueForChannel(chans[n]);

            if (requestedValue < maxVal)
            {
                fn->setCurrentChannel(chans[n]);
                fn->setParameter(FilterNode::PARAMETER_INDEX_LOW_CUT, requestedValue);
            }

            lastLowCutString = label->getText();
        } else if (label == orderValue) {
            if (requestedValue > 0 && roundDoubleToInt(requestedValue) <= FilterNode::MAX_ORDER) {
                fn->setCurrentChannel(chans[n]);
                fn->setParameter(FilterNode::PARAMETER_INDEX_ORDER, requestedValue);
                lastOrderString = label->getText();
            } else {
                // Revert any changes if it's an illegal value
                label->setText(lastOrderString, dontSendNotification);
            }
        }
    }

}

void FilterEditor::channelChanged (int channel, bool /*newState*/)
{
    FilterNode* fn = (FilterNode*) getProcessor();

    highCutValue->setText (String (fn->getHighCutValueForChannel (channel)), dontSendNotification);
    lowCutValue->setText  (String (fn->getLowCutValueForChannel  (channel)), dontSendNotification);
    orderValue->setText  (String (fn->getOrderValueForChannel(channel)), dontSendNotification);
    applyFilterOnChan->setToggleState (fn->getBypassStatusForChannel (channel), dontSendNotification);
}

void FilterEditor::buttonEvent(Button* button)
{

    if (button == applyFilterOnADC)
    {
        FilterNode* fn = (FilterNode*) getProcessor();
        fn->setApplyOnADC(applyFilterOnADC->getToggleState());

    }
    else if (button == applyFilterOnChan)
    {
        FilterNode* fn = (FilterNode*) getProcessor();

        Array<int> chans = getActiveChannels();

        for (int n = 0; n < chans.size(); n++)
        {
            float newValue = button->getToggleState() ? 1.0 : 0.0;

            fn->setCurrentChannel(chans[n]);
            fn->setParameter(FilterNode::PARAMETER_INDEX_ENABLE, newValue);
        }
    }
}


void FilterEditor::saveCustomParameters(XmlElement* xml)
{

    xml->setAttribute("Type", "FilterEditor");

    lastHighCutString = highCutValue->getText();
    lastLowCutString = lowCutValue->getText();
    lastOrderString = orderValue->getText();

    XmlElement* textLabelValues = xml->createNewChildElement("VALUES");
    textLabelValues->setAttribute("HighCut",lastHighCutString);
    textLabelValues->setAttribute("LowCut", lastLowCutString);
    textLabelValues->setAttribute("Order", lastOrderString);
    textLabelValues->setAttribute("ApplyToADC",	applyFilterOnADC->getToggleState());
}

void FilterEditor::loadCustomParameters(XmlElement* xml)
{

    forEachXmlChildElement(*xml, xmlNode)
    {
        if (xmlNode->hasTagName("VALUES"))
        {
            lastHighCutString = xmlNode->getStringAttribute("HighCut", lastHighCutString);
            lastLowCutString = xmlNode->getStringAttribute("LowCut", lastLowCutString);
            lastOrderString = xmlNode->getStringAttribute("Order", lastOrderString);
            resetToSavedText();

            applyFilterOnADC->setToggleState(xmlNode->getBoolAttribute("ApplyToADC",false), sendNotification);
        }
    }


}
