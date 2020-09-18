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

#include "HeartbeatNode.h"
#include "HeartbeatNodeEditor.h"
#include <memory>


HeartbeatNode::HeartbeatNode() : GenericProcessor("Heartbeat") {
    setProcessorType(PROCESSOR_TYPE_FILTER);
    buffer_.allocate(1, true);
}


AudioProcessorEditor *HeartbeatNode::createEditor() {
    editor = new HeartbeatNodeEditor(this, true);
    return editor;
}

void HeartbeatNode::createEventChannels() {
    auto *chan = new EventChannel(EventChannel::INT64_ARRAY,
                                  1,
                                  1,
                                  getSourceNode()->getSampleRate(getNodeId()),
                                  this);
    chan->setName("Heartbeat");
    chan->setDescription("Heartbeat events");
    chan->setIdentifier("heartbeat");
    eventChannelArray.add(chan);
}

void HeartbeatNode::process(AudioSampleBuffer &buffer) {
    if (buffer.getNumChannels() <= 0) {
        return;
    }

    int64 last_sample_received = getTimestamp(0) + buffer.getNumSamples();
    if (last_sample_received - last_heartbeat_at_ > heartbeat_every_samples_) {
        // emit an event with just the timestamp.
        MetaDataValueArray metadata;
        const EventChannel *chan = getEventChannel(getEventChannelIndex(0, getNodeId()));
        buffer_[0] = last_sample_received;
        BinaryEventPtr event = BinaryEvent::createBinaryEvent<juce::int64>(chan,
                                                                           last_sample_received,
                                                                           buffer_.getData(),
                                                                           sizeof(juce::int64),
                                                                           metadata);
        addEvent(chan, event, 0);
        last_heartbeat_at_ = last_sample_received;
    }
}

int HeartbeatNode::heartbeatEverySamples() const {
    return heartbeat_every_samples_;
}

void HeartbeatNode::setHeartbeatEverySamples(int period) {
    heartbeat_every_samples_ = period;
}

void HeartbeatNode::saveCustomParametersToXml(XmlElement *parentElement) {
    XmlElement *mainNode = parentElement->createNewChildElement("HeartbeatNode");
    mainNode->setAttribute("heartbeat_every_samples", heartbeatEverySamples());
}

void HeartbeatNode::loadCustomParametersFromXml() {
    if (parametersAsXml == nullptr) {
        return;
    }

    forEachXmlChildElement(*parametersAsXml, mainNode) {
        if (!mainNode->hasTagName("HeartbeatNode")) {
            continue;
        }

        heartbeat_every_samples_ = mainNode->getIntAttribute("heartbeat_every_samples", 1);
    }

    ((HeartbeatNodeEditor *) editor.get())->refreshFromProcessor();
}

