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

#include "HeartbeatNodeEditor.h"

HeartbeatNodeEditor::HeartbeatNodeEditor(
        GenericProcessor *parentNode, bool useDefaultParameterEditors)
        : GenericEditor(parentNode, useDefaultParameterEditors) {
    desiredWidth = 150;

    heartbeatEverySamplesLabel = newStaticLabel("Heartbeat Period", 10, 25, 80, 20);
    heartbeatEverySamplesLabelValue = newInputLabel("heartbeatEverySamplesLabelValue",
                                                    "Set the number of samples between each heartbeat",
                                                    15,
                                                    42,
                                                    80,
                                                    18);
    heartbeatEverySamplesLabelValue->addListener(this);
    refreshFromProcessor();
}

void HeartbeatNodeEditor::updateProcessor() {
    auto period = heartbeatEverySamplesLabelValue->getText().getIntValue();
    if (period > 0) {
        auto processor = dynamic_cast<HeartbeatNode *>(getProcessor());
        processor->setHeartbeatEverySamples(period);
    }
}


void HeartbeatNodeEditor::labelTextChanged(Label *label) {
    updateProcessor();
}

void HeartbeatNodeEditor::refreshFromProcessor() {
    auto processor = (HeartbeatNode *) getProcessor();
    heartbeatEverySamplesLabelValue->setText(juce::String(processor->heartbeatEverySamples()), dontSendNotification);
}
