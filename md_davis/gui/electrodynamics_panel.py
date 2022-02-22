import wx
from wx.lib.scrolledpanel import ScrolledPanel
from wx.lib.agw.floatspin import FloatSpin
import md_davis


class ElectrodynamicsPanel(wx.lib.scrolledpanel.ScrolledPanel):

    def __init__(self, parent):
        ScrolledPanel.__init__(self, parent)

        self.panel_sizer = wx.BoxSizer(orient=wx.VERTICAL)

        self.electrodynamics_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.panel_sizer.Add(self.electrodynamics_sizer, 0, wx.EXPAND, 5)
        self.electrodynamics_label = wx.StaticText(
            self, label='Electrostatics Directory')
        self.electrodynamics_sizer.Add(window=self.electrodynamics_label,
                                      proportion=0, flag=wx.ALL, border=5)
        self.electrodynamics_picker = wx.DirPickerCtrl(
            self,
            message='Select directory containing MD DaVis electrostatic files')
        self.electrodynamics_sizer.Add(window=self.electrodynamics_picker,
                                      proportion=1, flag=wx.ALL, border=5)

        self.name_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.panel_sizer.Add(self.name_sizer, 0, wx.EXPAND, 5)
        self.text_name = wx.StaticText(self, label='Name',)
        self.name_sizer.Add(window=self.text_name, proportion=0, flag=wx.ALL, border=5)
        self.name_textbox = wx.TextCtrl(self)
        self.name_sizer.Add(self.name_textbox, 1, wx.ALL, 5)

        self.time_step_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.panel_sizer.Add(self.time_step_sizer, 0, wx.EXPAND, 5)
        self.time_step_label = wx.StaticText(self, label='Time Step',
                                        )
        self.time_step = wx.SpinCtrl(self, min=1, max= 10000, initial=1)
        self.time_step_sizer.Add(self.time_step_label, 0, wx.ALL, 5)
        self.time_step_sizer.Add(self.time_step, 1, wx.ALL, 5)

        self.checkbox_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.panel_sizer.Add(self.checkbox_sizer, 0, wx.EXPAND, 5)
        self.pymol_checkbox = wx.CheckBox(self, label='Show PyMOL Window')
        self.surface_checkbox = wx.CheckBox(self, label='Show Surface')
        self.surface_checkbox.SetValue(True)
        self.secstr_checkbox = wx.CheckBox(self,
                                           label='Color by Secondary Structure')
        self.secstr_checkbox.SetValue(True)
        self.checkbox_sizer.Add(self.pymol_checkbox, 1, wx.ALL, 5)
        self.checkbox_sizer.Add(self.surface_checkbox, 1, wx.ALL, 5)
        self.checkbox_sizer.Add(self.secstr_checkbox, 1, wx.ALL, 5)

        self.line_space_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.panel_sizer.Add(self.line_space_sizer, 0, wx.EXPAND, 5)
        self.line_space_label = wx.StaticText(self, label='Spacing')
        self.line_spacing = FloatSpin(self, min_val=0, value=4, digits=2)
        self.line_length_label = wx.StaticText(self, label='Length')
        self.line_length = FloatSpin(self, min_val=0, value=10, digits=2)
        self.line_thickness_label = wx.StaticText(self, label='Thickness')
        self.line_thickness = FloatSpin(self, min_val=0, value=1, digits=2)
        self.line_space_sizer.Add(self.line_space_label, 0, wx.ALL, 5)
        self.line_space_sizer.Add(self.line_spacing, 1, wx.ALL, 5)
        self.line_space_sizer.Add(self.line_length_label, 0, wx.ALL, 5)
        self.line_space_sizer.Add(self.line_length, 1, wx.ALL, 5)
        self.line_space_sizer.Add(self.line_thickness_label, 0, wx.ALL, 5)
        self.line_space_sizer.Add(self.line_thickness, 1, wx.ALL, 5)

        self.theme_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.panel_sizer.Add(self.theme_sizer, 0, wx.EXPAND, 5)
        self.theme_label = wx.StaticText(self, label='Select Theme',
                                    )
        self.dark_radio = wx.RadioButton(self, label='Dark')
        self.dark_radio.SetValue(True)
        self.light_radio = wx.RadioButton(self, label='Light')
        self.theme_sizer.Add(self.theme_label, 0, wx.ALL, 5)
        self.theme_sizer.Add(self.dark_radio, 0, wx.ALL, 5)
        self.theme_sizer.Add(self.light_radio, 0, wx.ALL, 5)

        self.output_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.panel_sizer.Add(self.output_sizer, 0, wx.EXPAND, 5)
        self.output_label = wx.StaticText(self, label="Output File")
        self.output_sizer.Add(window=self.output_label, proportion=0, flag=wx.ALL, border=5)
        self.output_picker = wx.FilePickerCtrl(
            self, path='', message="Name for PyMOL session file",
            wildcard='PyMOL Session (*.pse)|*.pse',
            style=wx.FLP_SAVE | wx.FLP_CHANGE_DIR | wx.FLP_USE_TEXTCTRL)
        self.output_sizer.Add(window=self.output_picker, proportion=1, flag=wx.ALL, border=5)

        self.button_sizer = wx.BoxSizer(orient=wx.VERTICAL)
        self.button = wx.Button(self, label='Show Electric Field Dynamics')
        self.button.Bind(wx.EVT_BUTTON, self.on_button_press)
        self.button_sizer.Add(self.button, proportion=0,
                              flag=wx.ALL | wx.CENTER, border=5)
        self.panel_sizer.Add(self.button_sizer, 0, wx.EXPAND | wx.ALL, 5)

        self.SetSizer(self.panel_sizer)
        self.Layout()
        self.Fit()

    def on_button_press(self, event):
        """
        Animate the electric field dynamics in PyMOL
        @param event: The event object
        """

        if self.output_picker.GetPath() or self.pymol_checkbox.GetValue():
            md_davis.electrostatics.electrodynamics.get_electrodynamics(
                electrostatics_directory=self.electrodynamics_picker.GetPath(),
                name=self.name_textbox.GetValue(),
                surface=self.surface_checkbox.GetValue(),
                ss_color=self.secstr_checkbox.GetValue(),
                spacing=self.line_spacing.GetValue(),
                time_step=self.time_step.GetValue(),
                length=self.line_length.GetValue(),
                width=self.line_thickness.GetValue(),
                output=self.output_picker.GetPath(),
                hide=self.pymol_checkbox.GetValue(),
                light=self.light_radio.GetValue(),
            )
        else:
            wx.MessageBox("Please provide output file to save or select 'Show PyMOL Window'", "Message",
                          wx.OK | wx.ICON_INFORMATION)
