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

#ifndef __HEARTBEATNODEEDITOR_H_28EB4CC9__
#define __HEARTBEATNODEEDITOR_H_28EB4CC9__


#include <EditorHeaders.h>
#include <VisualizerEditorHeaders.h>
#include <VisualizerWindowHeaders.h>
#include "HeartbeatNode.h"

/**

  User interface for the HeartbeatNode processor.

  @see HeartbeatNode
*/
class HeartbeatNodeEditor : public GenericEditor,
                            public Label::Listener {
public:
    HeartbeatNodeEditor(GenericProcessor *parentNode, bool useDefaultParameterEditors);

    ~HeartbeatNodeEditor() override = default;

    // UI listeners
    void labelTextChanged(Label* label) override;

    // Non-override
    void refreshFromProcessor();

private:
    void updateProcessor();

    ScopedPointer<Label> heartbeatEverySamplesLabel;
    ScopedPointer<Label> heartbeatEverySamplesLabelValue;

    Label *newStaticLabel(
            const std::string& labelText,
            int boundsX,
            int boundsY,
            int boundsWidth,
            int boundsHeight) {
        return newStaticLabel(labelText, boundsX, boundsY, boundsWidth, boundsHeight, this);
    }

    static Label *newStaticLabel(
            const std::string& labelText,
            int boundsX,
            int boundsY,
            int boundsWidth,
            int boundsHeight,
            Component *parent) {
        auto *label = new Label(labelText, labelText);
        label->setBounds(boundsX, boundsY, boundsWidth, boundsHeight);
        label->setFont(Font("Small Text", 12, Font::plain));
        label->setColour(Label::textColourId, Colours::darkgrey);
        parent->addAndMakeVisible(label);
        return label;
    }

    Label *newInputLabel(
            const std::string &componentName,
            const std::string &tooltip,
            int boundsX,
            int boundsY,
            int boundsWidth,
            int boundsHeight) {
        return newInputLabel(componentName, tooltip, boundsX, boundsY, boundsWidth, boundsHeight, this);
    }

    static Label *newInputLabel(
            const std::string &componentName,
            const std::string &tooltip,
            int boundsX,
            int boundsY,
            int boundsWidth,
            int boundsHeight,
            Component *parent) {
        auto *label = new Label(componentName, "");
        label->setBounds(boundsX, boundsY, boundsWidth, boundsHeight);
        label->setFont(Font("Default", 15, Font::plain));
        label->setColour(Label::textColourId, Colours::white);
        label->setColour(Label::backgroundColourId, Colours::grey);
        label->setEditable(true);
        label->setTooltip(tooltip);
        parent->addAndMakeVisible(label);
        return label;
    }

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(HeartbeatNodeEditor)
};


#endif  // __HEARTBEATNODEEDITOR_H_28EB4CC9__
