import collections
import os
import wx
from wx.lib.scrolledpanel import ScrolledPanel
from wx.lib.intctrl import IntCtrl
from wx.lib.filebrowsebutton import FileBrowseButton
import md_davis


class LandscapePanel(wx.lib.scrolledpanel.ScrolledPanel):

    def __init__(self, parent):
        ScrolledPanel.__init__(self, parent)

        self.panel_sizer = wx.BoxSizer(orient=wx.VERTICAL)

        self.horrizontal_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.panel_sizer.Add(self.horrizontal_sizer, 0, wx.EXPAND, 5)

        self.box = wx.StaticBox(self, label='Data 1')
        self.box_sizer = wx.StaticBoxSizer(self.box, orient=wx.VERTICAL)
        self.horrizontal_sizer.Add(self.box_sizer, 0, wx.EXPAND, 5)

        self.X_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.X_label = wx.StaticText(self, label="X")
        self.X_sizer.Add(window=self.X_label, proportion=0, flag=wx.ALL, border=10)
        self.X_picker = wx.FilePickerCtrl(self, path='', message="Select a structure",
                                                  wildcard='Structure files (*.pdb;*.gro)|*.pdb;*.gro',
                                                  style=wx.FLP_CHANGE_DIR | wx.FLP_USE_TEXTCTRL)
        self.X_sizer.Add(window=self.X_picker, proportion=1, flag=wx.ALL, border=5)
        self.box_sizer.Add(self.X_sizer, 1, wx.EXPAND, 5)

        self.Y_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.Y_label = wx.StaticText(self, label="Y")
        self.Y_sizer.Add(window=self.Y_label, proportion=0, flag=wx.ALL, border=10)
        self.Y_picker = wx.FilePickerCtrl(self, path='', message="Select a structure",
                                                   wildcard='Structure files (*.trr;*.xtc)|*.trr;*.xtc',
                                                   style=wx.FLP_CHANGE_DIR | wx.FLP_USE_TEXTCTRL)
        self.Y_sizer.Add(window=self.Y_picker, proportion=1, flag=wx.ALL, border=5)
        self.box_sizer.Add(self.Y_sizer, 1, wx.EXPAND, 5)

        self.name_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.text_name = wx.StaticText(self, label="Name")
        self.name_sizer.Add(window=self.text_name, proportion=0, flag=wx.ALL, border=10)
        self.name_textbox = wx.TextCtrl(self)
        self.name_sizer.Add(self.name_textbox, 1, wx.ALL, 5)
        self.text_label = wx.StaticText(self, label="Label")
        self.name_sizer.Add(window=self.text_label, proportion=0, flag=wx.ALL, border=10)
        self.label_textbox = wx.TextCtrl(self)
        self.name_sizer.Add(self.label_textbox, 1, wx.ALL, 5)
        self.box_sizer.Add(self.name_sizer, 1, wx.EXPAND, 5)

        self.buttons_sizer = wx.BoxSizer(orient=wx.VERTICAL)
        self.horrizontal_sizer.Add(self.buttons_sizer, 1, wx.EXPAND, 5)
        self.buttons_sizer.Add((0, 0), proportion=1)
        self.del_button = wx.Button(self, label='-', size=(25, -1))
        self.buttons_sizer.Add(self.del_button, flag=wx.ALL, border=5)
        self.add_button = wx.Button(self, label='+', size=(25, -1))
        self.buttons_sizer.Add(self.add_button, flag=wx.ALL, border=5)

        self.SetSizer(self.panel_sizer)
        self.SetupScrolling()
        self.Layout()

    def onAddWidget(self, event):
        """"""
        self.number_of_buttons += 1
        label = "Button %s" %  self.number_of_buttons
        name = "button%s" % self.number_of_buttons
        new_button = wx.Button(self, label=label, name=name)
        self.widgetSizer.Add(new_button, 0, wx.ALL, 5)
        self.frame.fSizer.Layout()
        self.frame.Fit()

    def onRemoveWidget(self, event):
        if self.widgetSizer.GetChildren():
            sizer_item = self.widgetSizer.GetItem(self.number_of_buttons-1)
            widget = sizer_item.GetWindow()
            self.widgetSizer.Hide(widget)
            widget.Destroy()
            self.number_of_buttons -= 1
            self.frame.fSizer.Layout()
            self.frame.Fit()


    # def on_choice(self, event):
    #     """
    #     Add label box on selecting FASTA format
    #     """
    #     if self.dropdown.GetSelection() == 1:
    #         self.label_for_fasta.Show()
    #         self.label_box.Show()
    #         self.Layout()
    #     else:
    #         self.label_for_fasta.Show(False)
    #         self.label_box.Show(False)
    #         self.Layout()
    #
    # def on_button_press(self, event):
    #     """
    #     Browse for a PDB file
    #     @param event: The event object
    #     """
    #     output = ['toml', 'fasta', 'modeller', 'dict']
    #     pdb_file = self.file_picker.GetPath()
    #     if pdb_file:
    #         seq = md_davis.sequence.get_sequence(
    #             structure=pdb_file,
    #             label=self.label_box.GetValue(),
    #             return_type=output[self.dropdown.GetSelection()]
    #         )
    #         self.output_box.SetValue(str(seq))
