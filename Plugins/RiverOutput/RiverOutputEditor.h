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
#include "RiverOutput.h"

class ImageIcon;

/**

  User interface for the RiverOutput processor.

  @see RiverOutput
*/

class RiverOutputEditor : public GenericEditor,
                          public Label::Listener {
public:
    RiverOutputEditor(GenericProcessor *parentNode, bool useDefaultParameterEditors);

    ~RiverOutputEditor() override = default;

    void labelTextChanged(Label* label) override;

    void refreshLabels();

private:
    ImageIcon* icon;

    Label *newStaticLabel(
            const std::string& labelText,
            int boundsX,
            int boundsY,
            int boundsWidth,
            int boundsHeight);

    Label *newInputLabel(
            const std::string& componentName,
            const std::string& tooltip,
            int boundsX,
            int boundsY,
            int boundsWidth,
            int boundsHeight);

    ScopedPointer<Label> streamNameLabel;
    ScopedPointer<Label> streamNameLabelValue;

    ScopedPointer<Label> hostnameLabel;
    ScopedPointer<Label> hostnameLabelValue;

    ScopedPointer<Label> portLabel;
    ScopedPointer<Label> portLabelValue;
    std::string lastPortValue;

    ScopedPointer<Label> passwordLabel;
    ScopedPointer<Label> passwordLabelValue;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(RiverOutputEditor)
};




#endif  // __RIVEROUTPUTEDITOR_H_28EB4CC9__
