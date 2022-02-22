import wx
from wx.lib.scrolledpanel import ScrolledPanel
from wx.lib.agw.floatspin import FloatSpin
import md_davis


class LandscapePanel(wx.lib.scrolledpanel.ScrolledPanel):

    def __init__(self, parent):
        ScrolledPanel.__init__(self, parent)

        self.x_picker_array = []
        self.y_picker_array = []
        self.name_array = []
        self.label_array = []

        self.panel_sizer = wx.BoxSizer(orient=wx.VERTICAL)

        self.panel_sizer.Add(self.top_sizer(), 0, wx.EXPAND | wx.ALL, 5)
        self.panel_sizer.Add(self.input_sizer(), 0, wx.EXPAND, 5)

        self.SetSizer(self.panel_sizer)
        self.SetupScrolling()
        self.Layout()

    def input_sizer(self):
        self.horizontal_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)

        self.box = wx.StaticBox(self, label='Input')
        self.box_sizer = wx.StaticBoxSizer(self.box, orient=wx.VERTICAL)
        self.horizontal_sizer.Add(self.box_sizer, 1, wx.EXPAND, 5)

        self.X_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.X_label = wx.StaticText(self, label="X")
        self.X_sizer.Add(window=self.X_label, proportion=0,
                         flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL,
                         border=5)
        self.x_picker = wx.FilePickerCtrl(
            self, path='',
            message="Select a structure",
            style=wx.FLP_CHANGE_DIR | wx.FLP_USE_TEXTCTRL)
        self.X_sizer.Add(window=self.x_picker, proportion=1, flag=wx.ALL,
                         border=5)
        self.box_sizer.Add(self.X_sizer, 1, wx.EXPAND, 5)

        self.Y_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.Y_label = wx.StaticText(self, label="Y")
        self.Y_sizer.Add(window=self.Y_label, proportion=0,
                         flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL, border=5)
        self.y_picker = wx.FilePickerCtrl(
            self, path='',
            message="Select a structure",
            style=wx.FLP_CHANGE_DIR | wx.FLP_USE_TEXTCTRL)
        self.Y_sizer.Add(window=self.y_picker, proportion=1, flag=wx.ALL,
                         border=5)
        self.box_sizer.Add(self.Y_sizer, 1, wx.EXPAND, 5)

        self.name_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.text_name = wx.StaticText(self, label="Name")
        self.name_sizer.Add(window=self.text_name, proportion=0,
                            flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL,
                            border=5)
        self.name_textbox = wx.TextCtrl(self)
        self.name_sizer.Add(self.name_textbox, 1, wx.ALL, 5)
        self.text_label = wx.StaticText(self, label="Label")
        self.name_sizer.Add(window=self.text_label, proportion=0,
                            flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL,
                            border=5)
        self.label_textbox = wx.TextCtrl(self)
        self.name_sizer.Add(self.label_textbox, 1, wx.ALL, 5)
        self.box_sizer.Add(self.name_sizer, 1, wx.EXPAND, 5)

        self.buttons_sizer = wx.BoxSizer(orient=wx.VERTICAL)
        self.horizontal_sizer.Add(self.buttons_sizer, 0, wx.EXPAND, 5)
        self.buttons_sizer.Add((0, 0), proportion=1)
        self.del_button = wx.Button(self, label='-', size=(25, -1))
        self.del_button.Bind(wx.EVT_BUTTON, self.on_del_button)
        self.buttons_sizer.Add(self.del_button, flag=wx.ALL, border=5)
        self.add_button = wx.Button(self, label='+', size=(25, -1))
        self.add_button.Bind(wx.EVT_BUTTON, self.on_add_button)
        self.buttons_sizer.Add(self.add_button, flag=wx.ALL, border=5)

        self.x_picker_array.append(self.x_picker)
        self.y_picker_array.append(self.y_picker)
        self.name_array.append(self.name_textbox)
        self.label_array.append(self.label_textbox)
        return self.horizontal_sizer

    def top_sizer(self):
        self.run_sizer = wx.BoxSizer(orient=wx.VERTICAL)

        self.type_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.run_sizer.Add(self.type_sizer, 0, wx.ALL, 5)
        self.hist_radio = wx.RadioButton(self, label='Histogram',
                                         style=wx.RB_GROUP)
        self.hist_radio.Bind(wx.EVT_RADIOBUTTON, self.on_radio_choice)
        self.type_sizer.Add(self.hist_radio, 0,
                            wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
        self.fel_radio = wx.RadioButton(self, label='Free Energy')
        self.fel_radio.Bind(wx.EVT_RADIOBUTTON, self.on_radio_choice)
        self.fel_radio.SetValue(True)
        self.type_sizer.Add(self.fel_radio, 0,
                            wx.TOP | wx.BOTTOM | wx.ALIGN_CENTER_VERTICAL, 5)
        self.temperature_label = wx.StaticText(self, label='Temperature (in K)')
        self.temperature = FloatSpin(self, min_val=0, value=300, digits=2)
        self.type_sizer.Add(self.temperature_label, 0,
                            wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.type_sizer.Add(self.temperature, 0, wx.ALL, 5)

        self.view_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.run_sizer.Add(self.view_sizer, 0, wx.EXPAND, 5)
        self.perspective_radio = wx.RadioButton(self, label='Perspective',
                                                style=wx.RB_GROUP)
        self.perspective_radio.SetValue(True)
        self.view_sizer.Add(self.perspective_radio, 0,
                            wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.orthographic_radio = wx.RadioButton(self, label='Orthographic')
        self.view_sizer.Add(self.orthographic_radio, 0,
                            wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)

        self.title_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.run_sizer.Add(self.title_sizer, proportion=1,
                           flag=wx.EXPAND | wx.ALL,
                           border=5)
        self.title_text = wx.StaticText(self, label="Title")
        self.title_sizer.Add(window=self.title_text, proportion=0,
                             flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL,
                             border=5)
        self.title_textbox = wx.TextCtrl(self)
        self.title_sizer.Add(self.title_textbox, 1,
                             flag=wx.ALL | wx.ALIGN_CENTER_VERTICAL,
                             border=5)

        self.bins_box = wx.StaticBox(self, label='Shape of each landscape')
        self.bins_sizer = wx.StaticBoxSizer(self.bins_box, orient=wx.HORIZONTAL)
        self.run_sizer.Add(self.bins_sizer, proportion=1,
                           flag=wx.EXPAND | wx.ALL,
                           border=5)
        self.x_bins_text = wx.StaticText(self, label="X-bins")
        self.bins_sizer.Add(window=self.x_bins_text, proportion=0,
                            flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL,
                            border=5)
        self.x_bins = wx.SpinCtrl(self, min=1, max=1000, initial=100,
                                  size=(85, -1))
        self.bins_sizer.Add(self.x_bins, proportion=1,
                            flag=wx.ALL,
                            border=5)
        self.y_bins_text = wx.StaticText(self, label="Y-bins")
        self.bins_sizer.Add(window=self.y_bins_text, proportion=0,
                            flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL,
                            border=20)
        self.y_bins = wx.SpinCtrl(self, min=1, max=1000, initial=100,
                                  size=(85, -1))
        self.bins_sizer.Add(self.y_bins, proportion=1,
                            flag=wx.ALL,
                            border=5)
        self.common_checkbox = wx.CheckBox(self, label='Use Common Ranges')
        self.common_checkbox.SetValue(True)
        self.bins_sizer.Add(self.common_checkbox, proportion=0,
                            flag=wx.LEFT | wx.ALIGN_CENTER_VERTICAL,
                            border=20)

        self.axis_box = wx.StaticBox(self, label='Axis Labels')
        self.axis_sizer = wx.StaticBoxSizer(self.axis_box, orient=wx.HORIZONTAL)
        self.run_sizer.Add(self.axis_sizer, proportion=1,
                           flag=wx.EXPAND | wx.ALL, border=5)
        self.x_text = wx.StaticText(self, label="X")
        self.axis_sizer.Add(window=self.x_text, proportion=0, flag=wx.ALL,
                            border=10)
        self.x_label = wx.TextCtrl(self, value="RMSD (in nm)")
        self.axis_sizer.Add(self.x_label, 1, wx.ALL, 5)
        self.y_text = wx.StaticText(self, label="Y")
        self.axis_sizer.Add(window=self.y_text, proportion=0, flag=wx.ALL,
                            border=10)
        self.y_label = wx.TextCtrl(self, value="Radius of Gyration (in nm)")
        self.axis_sizer.Add(self.y_label, 1, wx.ALL, 5)
        self.z_text = wx.StaticText(self, label="Z")
        self.axis_sizer.Add(window=self.z_text, proportion=0, flag=wx.ALL,
                            border=10)
        self.z_label = wx.TextCtrl(self, value="Free Envergy (in kJ mol<sup>-1</sup>)")
        self.axis_sizer.Add(self.z_label, 1, wx.ALL, 5)

        self.plot_shape_box = wx.StaticBox(self,
                                           label='Size and layout of the plot')
        self.plot_shape_sizer = wx.StaticBoxSizer(self.plot_shape_box,
                                                  orient=wx.HORIZONTAL)
        self.run_sizer.Add(self.plot_shape_sizer, proportion=1,
                           flag=wx.EXPAND | wx.ALL, border=5)

        self.width_text = wx.StaticText(self, label="Width")
        self.plot_shape_sizer.Add(window=self.width_text, proportion=0,
                                  flag=wx.ALL, border=10)
        self.width = wx.SpinCtrl(self, value='', min=800, max=10000, initial=1920)
        self.plot_shape_sizer.Add(self.width, 1, wx.ALL, 5)
        self.height_text = wx.StaticText(self, label="Height")
        self.plot_shape_sizer.Add(window=self.height_text, proportion=0,
                                  flag=wx.ALL, border=10)
        self.height = wx.SpinCtrl(self, value='', min=600, max=10000, initial=1080)
        self.plot_shape_sizer.Add(self.height, 1, wx.ALL, 5)

        self.save_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.run_sizer.Add(self.save_sizer, 0, wx.EXPAND, 5)
        self.save_label = wx.StaticText(self, label='Save landscapes to')
        self.save_sizer.Add(self.save_label, 0,
                            wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.save_picker = wx.FilePickerCtrl(
            self, path='', message="Select HDF file to save the landscapes",
            wildcard='HDF file (*.hdf;*.h5;*.hdf5;*.he5)|*.h5;*.hdf5;*.he5;*.hdf',
            style=wx.FLP_SAVE | wx.FLP_CHANGE_DIR | wx.FLP_USE_TEXTCTRL | wx.FLP_OVERWRITE_PROMPT)
        self.save_sizer.Add(self.save_picker, 1, wx.ALL, 5)

        self.output_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.run_sizer.Add(self.output_sizer, 0, wx.EXPAND, 5)
        self.output_label = wx.StaticText(self, label='Output HTML file')
        self.output_sizer.Add(self.output_label, 0,
                              wx.ALL | wx.ALIGN_CENTER_VERTICAL, 5)
        self.output_picker = wx.FilePickerCtrl(
            self, path='', message="Select name for output HTML file",
            wildcard='HTML file (*.html)|*.html',
            style=wx.FLP_SAVE | wx.FLP_CHANGE_DIR | wx.FLP_USE_TEXTCTRL | wx.FLP_OVERWRITE_PROMPT)
        self.output_sizer.Add(self.output_picker, 1, wx.ALL, 5)

        self.button = wx.Button(self, label='Create Free Energy Landscape')
        self.button.Bind(wx.EVT_BUTTON, self.on_button_press)
        self.run_sizer.Add(self.button, proportion=0,
                           flag=wx.ALL | wx.CENTER, border=5)
        return self.run_sizer

    def on_radio_choice(self, event):
        """ Show temperature box if free energy landscape is selected """
        if self.fel_radio.GetValue():
            self.temperature_label.Show()
            self.temperature.Show()
            self.Layout()
        else:
            self.temperature_label.Show(False)
            self.temperature.Show(False)
            self.Layout()

    def on_add_button(self, event):
        """ Add a new input panel at the bottom """

        self.panel_sizer.Add(self.input_sizer(), 0, wx.EXPAND, 5)
        self.SetSizer(self.panel_sizer)
        self.SetupScrolling(scrollToTop=False)
        self.Layout()

    def on_del_button(self, event):
        """ Delete the last input panel """
        num_items = self.panel_sizer.GetItemCount()
        if num_items > 2:
            second_last_sizer = self.panel_sizer.GetChildren()[-1]
            second_last_sizer.DeleteWindows()
            self.panel_sizer.Remove(num_items - 1)
            self.SetSizer(self.panel_sizer)
            self.SetupScrolling(scrollToTop=False)
            self.Layout()
            self.x_picker_array.pop()
            self.y_picker_array.pop()
            self.name_array.pop()
            self.label_array.pop()

    def on_button_press(self, event):
        """
        Browse for a PDB file
        @param event: The event object
        """

        if self.save_picker.GetPath():
            save = self.save_picker.GetPath()
        else:
            save = None

        if self.output_picker.GetPath():
            output = self.output_picker.GetPath()
        else:
            wx.MessageBox("Please provide a name for the output HTML file",
                          "Message", wx.OK | wx.ICON_INFORMATION)
            return

        if self.fel_radio.GetValue():
            temperature = self.temperature.GetValue()
        else:
            temperature = None

        md_davis.landscape.landscape_xvg.landscape_xvg(
            x=[_.GetPath() for _ in self.x_picker_array],
            y=[_.GetPath() for _ in self.y_picker_array],
            name=[_.GetValue() for _ in self.name_array],
            label=[_.GetValue() for _ in self.label_array],
            temperature=temperature,
            common=self.common_checkbox.GetValue(),
            output=output,
            save=save,
            title=self.title_textbox.GetValue(),
            shape=(self.x_bins.GetValue(), self.y_bins.GetValue()),
            begin=0,
            end=None,
            limits=None,
            orthographic=self.orthographic_radio.GetValue(),
            layout=None,
            width=self.width.GetValue(),
            height=self.height.GetValue(),
            font=None,
            font_size=None,
            dtick=None,
            axis_labels="{ 'x': '" + self.x_label.GetValue() +
                        "', 'y': '" + self.y_label.GetValue() +
                        "', 'z': '" + self.z_label.GetValue() + "' }",
            columns=None
        )
        if self.hist_radio.GetValue():
            print(f'Histograms plotted in {output}')
        else:
            print(f'Free energy landscape plotted in {output}')
