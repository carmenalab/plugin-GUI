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

#include "RiverOutput.h"
#include "RiverOutputEditor.h"
#include <nlohmann/json.hpp>
#include <memory>

using json = nlohmann::json;


RiverOutput::RiverOutput()
        : GenericProcessor("River Output"),
          spike_schema_({river::FieldDefinition("channel_index", river::FieldDefinition::INT32, 4),
                         river::FieldDefinition("unit_index", river::FieldDefinition::INT32, 4),
                         river::FieldDefinition("data_index", river::FieldDefinition::INT64, 8),
                        }) {
    setProcessorType(PROCESSOR_TYPE_SINK);

    // Start with some sane defaults.
    redis_connection_hostname_ = "127.0.0.1";
    redis_connection_port_ = 6379;

    stream_name_ = new Parameter("stream_name", juce::String(""), 0, true);
    stream_name_->addListener(this);
    parameters.add(stream_name_);
}


AudioProcessorEditor *RiverOutput::createEditor() {
    editor = new RiverOutputEditor(this, true);
    return editor;
}


void RiverOutput::handleSpike(const SpikeChannel *spikeInfo, const MidiMessage &event, int) {
    SpikeEventPtr spike = SpikeEvent::deserializeFromMessage(event, spikeInfo);
    if (!spike) {
        return;
    }

    RiverSpike river_spike;
    river_spike.channel_index = getSpikeChannelIndex(spike);
    river_spike.data_index = spike->getTimestamp();
    river_spike.unit_index = ((int) spike->getSortedID()) - 1;
    writer_->Write(&river_spike, 1);
}

void RiverOutput::handleEvent(const EventChannel* eventInfo, const MidiMessage& msg, int) {
    BinaryEventPtr event = BinaryEvent::deserializeFromMessage(msg, eventInfo);
    if (!event) {
        return;
    }

    const char *ptr = (const char *) event->getBinaryDataPointer();
    size_t data_size = eventInfo->getDataSize();

    // Assert (when compiled in debug) that the sizes match up.
    jassert(((int) data_size == writer_->schema().sample_size()));

    // Assume that the binary data in the event matches the sample size exactly. If it doesn't, crashes will happen!
    writer_->WriteBytes(ptr, 1);
}

void RiverOutput::parameterValueChanged(Value &, const String &) {
    // must be stream name; we already use the Parameter as source of truth for the rest so just ping the editor
    if (editor) {
        const MessageManagerLock mm;
        ((RiverOutputEditor *) editor.get())->refreshLabelsFromProcessor();
    }
}

bool RiverOutput::enable() {
    auto sn = streamName();
    if (sn.empty() || redis_connection_hostname_.empty() || redis_connection_port_ <= 0) {
        return false;
    }

    jassert(!writer_);

    river::RedisConnection connection(
            redis_connection_hostname_,
            redis_connection_port_,
            redis_connection_password_);

    // TODO: what happens if invalid redis connection?
    auto new_writer = std::make_unique<river::StreamWriter>(connection);
    writer_.swap(new_writer);

    std::unordered_map<std::string, std::string> metadata;
    if (shouldConsumeSpikes()) {
        if (getTotalSpikeChannels() == 0) {
            // Can't consume spikes if there are no spike channels.
            return false;
        }

        // Assume that all spike channels have the same details.
        auto spike_channel = getSpikeChannel(0);
        metadata["prepeak_samples"] = std::to_string(spike_channel->getPrePeakSamples());
        metadata["postpeak_samples"] = std::to_string(spike_channel->getPostPeakSamples());
    }

    writer_->Initialize(sn, getSchema(), metadata);

    if (editor) {
        // GenericEditor#enable isn't marked as virtual, so need to *upcast* to VisualizerEditor :(
        ((VisualizerEditor *) (editor.get()))->enable();
    }

    return true;
}


bool RiverOutput::disable() {
    if (writer_) {
        writer_->Stop();
        writer_.reset();
    }
    return true;
}


void RiverOutput::process(AudioSampleBuffer &buffer) {
    if (writer_) {
        checkForEvents(shouldConsumeSpikes());
    }
}

std::string RiverOutput::streamName() const {
    return stream_name_->getValue(0).toString().toStdString();
}

int64_t RiverOutput::totalSamplesWritten() const {
    if (writer_) {
        return writer_->total_samples_written();
    } else {
        return 0;
    }
}

const std::string &RiverOutput::redisConnectionHostname() const {
    return redis_connection_hostname_;
}

void RiverOutput::setRedisConnectionHostname(const std::string &redisConnectionHostname) {
    redis_connection_hostname_ = redisConnectionHostname;
}

int RiverOutput::redisConnectionPort() const {
    return redis_connection_port_;
}

void RiverOutput::setRedisConnectionPort(int redisConnectionPort) {
    redis_connection_port_ = redisConnectionPort;
}

const std::string &RiverOutput::redisConnectionPassword() const {
    return redis_connection_password_;
}

void RiverOutput::setRedisConnectionPassword(const std::string &redisConnectionPassword) {
    redis_connection_password_ = redisConnectionPassword;
}

void RiverOutput::saveCustomParametersToXml(XmlElement *parentElement) {
    XmlElement *mainNode = parentElement->createNewChildElement("RiverOutput");
    mainNode->setAttribute("hostname", redisConnectionHostname());
    mainNode->setAttribute("port", redisConnectionPort());
    mainNode->setAttribute("password", redisConnectionPassword());

    if (event_schema_) {
        std::string event_schema_json = event_schema_->ToJson();
        mainNode->setAttribute("event_schema_json", event_schema_json);
    }
}

void RiverOutput::loadCustomParametersFromXml() {
    if (parametersAsXml == nullptr) {
        return;
    }

    forEachXmlChildElement(*parametersAsXml, mainNode) {
        if (!mainNode->hasTagName("RiverOutput")) {
            continue;
        }

        redis_connection_hostname_ = mainNode->getStringAttribute("hostname", "127.0.0.1").toStdString();
        redis_connection_port_ = mainNode->getIntAttribute("port", 6379);
        redis_connection_password_ = mainNode->getStringAttribute("password", "").toStdString();
        if (mainNode->hasAttribute("event_schema_json")) {
            String j = mainNode->getStringAttribute("event_schema_json");
            const river::StreamSchema& schema = river::StreamSchema::FromJson(j.toStdString());
            setEventSchema(schema);
        } else {
            clearEventSchema();
        }
    }

    ((RiverOutputEditor *) editor.get())->refreshSchemaFromProcessor();
    ((RiverOutputEditor *) editor.get())->refreshLabelsFromProcessor();
}

void RiverOutput::setEventSchema(const river::StreamSchema& eventSchema) {
    auto p = std::make_shared<river::StreamSchema>(eventSchema);
    event_schema_.swap(p);
}

void RiverOutput::clearEventSchema() {
    event_schema_.reset();
}

bool RiverOutput::shouldConsumeSpikes() const {
    return !event_schema_;
}

river::StreamSchema RiverOutput::getSchema() const {
    if (event_schema_) {
        return *event_schema_;
    }
    return spike_schema_;
}
