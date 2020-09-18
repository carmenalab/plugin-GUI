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

#ifndef __RIVEROUTPUTEDITOR_H_28EB4CC9__
#define __RIVEROUTPUTEDITOR_H_28EB4CC9__


#include <EditorHeaders.h>
#include <VisualizerEditorHeaders.h>
#include <VisualizerWindowHeaders.h>
#include "RiverOutput.h"
#include "SchemaListBox.h"

class ImageIcon;

class RiverOutputCanvas;


/**

  User interface for the RiverOutput processor.

  @see RiverOutput
*/
class RiverOutputEditor : public VisualizerEditor,
                          public Label::Listener,
                          public ComboBox::Listener {
public:
    RiverOutputEditor(GenericProcessor *parentNode, bool useDefaultParameterEditors);

    ~RiverOutputEditor() override = default;

    Visualizer* createNewCanvas() override;

    // UI listeners
    void labelTextChanged(Label* label) override;
    void comboBoxChanged(ComboBox *comboBoxThatHasChanged) override;
    void buttonEvent(Button* button) override;

    // Non-overrides:
    void refreshLabelsFromProcessor();
    void refreshSchemaFromProcessor();

    Component *getOptionsPanel() {
        return optionsPanel.get();
    }

private:
    ImageIcon* icon;

    void updateProcessorSchema();

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
            Component *parent);

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
            Component *parent);

    ScopedPointer<Label> hostnameLabel;
    ScopedPointer<Label> hostnameLabelValue;

    ScopedPointer<Label> portLabel;
    ScopedPointer<Label> portLabelValue;
    std::string lastPortValue;

    ScopedPointer<Label> passwordLabel;
    ScopedPointer<Label> passwordLabelValue;

    // OPTIONS PANEL
    RiverOutputCanvas* canvas;
    ScopedPointer<Component> optionsPanel;
    ScopedPointer<Label> optionsPanelTitle;
    ScopedPointer<Label> inputTypeTitle;

    ScopedPointer<Label> streamNameLabel;
    ScopedPointer<Label> streamNameLabelValue;

    ScopedPointer<Label> totalSamplesWrittenLabel;
    ScopedPointer<Label> totalSamplesWrittenLabelValue;

    // OPTIONS PANEL: Input Type
    const int inputTypeRadioId = 1;
    ScopedPointer<ToggleButton> inputTypeSpikeButton;
    ScopedPointer<ToggleButton> inputTypeEventButton;

    ScopedPointer<Label> fieldNameLabel;
    ScopedPointer<Label> fieldNameLabelValue;

    ScopedPointer<Label> fieldTypeLabel;
    ScopedPointer<ComboBox> fieldTypeComboBox;

    ScopedPointer<UtilityButton> addFieldButton;
    ScopedPointer<UtilityButton> removeSelectedFieldButton;

    ScopedPointer<SchemaListBox> schemaList;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(RiverOutputEditor)
};

class RiverOutputCanvas : public Visualizer {
public:
    explicit RiverOutputCanvas(GenericProcessor *n);
    ~RiverOutputCanvas() override = default;

    void refreshState() override;
    void update() override;
    void refresh() override;
    void beginAnimation() override;
    void endAnimation() override;
    void setParameter(int, float) override;
    void setParameter(int, int, int, float) override;
    void paint(Graphics& g) override;
    void resized() override;

private:
    ScopedPointer<Viewport> viewport;
    GenericProcessor* processor;
    RiverOutputEditor* editor;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(RiverOutputCanvas)
};


#endif  // __RIVEROUTPUTEDITOR_H_28EB4CC9__
