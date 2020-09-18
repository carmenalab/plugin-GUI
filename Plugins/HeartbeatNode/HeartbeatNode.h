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

#ifndef __HEARTBEATNODE_H_F7BDA585__
#define __HEARTBEATNODE_H_F7BDA585__

#include <ProcessorHeaders.h>


/**
 *  TODO
    @see GenericProcessor
 */
class HeartbeatNode : public GenericProcessor {
public:
    HeartbeatNode();

    ~HeartbeatNode() override = default;

    void createEventChannels() override;

    /** Searches for events and triggers the River output when appropriate. */
    void process(AudioSampleBuffer &buffer) override;

    /** Creates the HeartbeatNodeEditor. */
    AudioProcessorEditor *createEditor() override;

    /** To/from conversion XML */
    void saveCustomParametersToXml (XmlElement* parentElement) override;
    /** Load custom settings from XML*/
    void loadCustomParametersFromXml() override;

    //
    // Non-override methods:
    //
    int heartbeatEverySamples() const;
    void setHeartbeatEverySamples(int period);

private:
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (HeartbeatNode)

    int heartbeat_every_samples_ = 1;
    int64_t last_heartbeat_at_ = -1;
    HeapBlock<juce::int64> buffer_;
};


#endif  // __HEARTBEATNODE_H_F7BDA585__
