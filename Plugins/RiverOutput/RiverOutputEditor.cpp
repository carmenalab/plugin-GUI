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

#include "RiverOutputEditor.h"

RiverOutputEditor::RiverOutputEditor(
        GenericProcessor *parentNode, bool useDefaultParameterEditors)
        : GenericEditor(parentNode, useDefaultParameterEditors) {

    desiredWidth = 350;

    hostnameLabel = newStaticLabel("Hostname", 10,25,80,20);
    hostnameLabelValue = newInputLabel("hostnameLabelValue", "Set the hostname for River", 15,42,80,18);

    portLabel = newStaticLabel("Port", 10, 65, 80, 20);
    portLabelValue = newInputLabel("hostnamePortValue", "Set the port for River", 15, 82, 60, 18);

    passwordLabel = newStaticLabel("Password", 105, 25, 100, 20);
    passwordLabelValue = newInputLabel("hostnamePasswordValue", "Set the password for River", 110, 42, 100, 18);

    streamNameLabel = newStaticLabel("Stream Name", 105, 65, 80, 20);
    streamNameLabelValue = newInputLabel("streamNameLabelValue", "Set the port for River", 110, 82, 120, 18);
    streamNameLabelValue->setEditable(false, false);
    refreshLabels();
}

void RiverOutputEditor::labelTextChanged(Label *label) {
    auto river = (RiverOutput *) getProcessor();
    if (label == hostnameLabelValue) {
        river->setRedisConnectionHostname(label->getText().toStdString());
    } else if (label == portLabelValue) {
        int port = label->getText().getIntValue();
        if (port > 0) {
            river->setRedisConnectionPort(port);
            lastPortValue = label->getText().toStdString();
        } else {
            label->setText(juce::String(lastPortValue), dontSendNotification);
        }
    } else if (label == passwordLabelValue) {
        river->setRedisConnectionPassword(label->getText().toStdString());
    }
}

void RiverOutputEditor::refreshLabels() {
    auto river = (RiverOutput *) getProcessor();
    lastPortValue = std::to_string(river->redisConnectionPort());

    hostnameLabelValue->setText(river->redisConnectionHostname(), dontSendNotification);
    portLabelValue->setText(lastPortValue, dontSendNotification);
    passwordLabelValue->setText(river->redisConnectionPassword(), dontSendNotification);
    streamNameLabelValue->setText(river->streamName(), dontSendNotification);
}

Label *RiverOutputEditor::newStaticLabel(
        const std::string& labelText, int boundsX, int boundsY, int boundsWidth, int boundsHeight) {
    auto *label = new Label(labelText, labelText);
    label->setBounds(boundsX, boundsY, boundsWidth, boundsHeight);
    label->setFont(Font("Small Text", 12, Font::plain));
    label->setColour(Label::textColourId, Colours::darkgrey);
    addAndMakeVisible(label);
    return label;
}

Label *RiverOutputEditor::newInputLabel(
        const std::string& componentName,
        const std::string& tooltip,
        int boundsX,
        int boundsY,
        int boundsWidth,
        int boundsHeight) {
    auto *label = new Label(componentName, "");
    label->setBounds(boundsX, boundsY, boundsWidth, boundsHeight);
    label->setFont(Font("Default", 15, Font::plain));
    label->setColour(Label::textColourId, Colours::white);
    label->setColour(Label::backgroundColourId, Colours::grey);
    label->setEditable(true);
    label->addListener(this);
    label->setTooltip(tooltip);
    addAndMakeVisible(label);
    return label;
}
