import collections
import os
import wx
from wx.lib.scrolledpanel import ScrolledPanel
from wx.lib.intctrl import IntCtrl
from wx.lib.filebrowsebutton import FileBrowseButton
import md_davis


class ElectrostaticsPanel(wx.lib.scrolledpanel.ScrolledPanel):

    def __init__(self, parent):
        ScrolledPanel.__init__(self, parent)

        self.panel_sizer = wx.BoxSizer(orient=wx.VERTICAL)

        self.file_picker_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.input_label = wx.StaticText(self, label="PDB File")
        self.file_picker_sizer.Add(window=self.input_label, proportion=0, flag=wx.ALL, border=10)
        self.wildcard = 'PDB files (*.pdb)|*.pdb'
        self.file_picker = wx.FilePickerCtrl(self, path='',
                                             message="Select a PDB file",
                                             wildcard=self.wildcard)
        self.file_picker_sizer.Add(window=self.file_picker, proportion=1, flag=wx.ALL, border=5)

        self.run_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.select_output = wx.StaticText(self, label="Select Output", size=(85, -1))
        self.run_sizer.Add(window=self.select_output, proportion=0, flag=wx.ALL, border=10)

        self.output_type = ['TOML', 'FASTA', 'MODELLER', 'Python Dictionary']
        self.dropdown = wx.Choice(self, choices=self.output_type)
        self.dropdown.SetSelection(0)
        self.run_sizer.Add(self.dropdown, proportion=1, flag=wx.ALL | wx.CENTER, border=5)
        self.dropdown.Bind(wx.EVT_CHOICE, self.on_choice)

        self.button = wx.Button(self, label='Get Sequence')
        self.button.Bind(wx.EVT_BUTTON, self.on_button_press)
        self.run_sizer.Add(self.button, proportion=0, flag=wx.ALL | wx.CENTER, border=5)

        self.label_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.label_for_fasta = wx.StaticText(self, label="Label for FASTA", size=(85, -1))
        self.label_sizer.Add(window=self.label_for_fasta, proportion=0, flag=wx.ALL, border=10)
        self.label_box = wx.TextCtrl(self)
        self.label_sizer.Add(self.label_box, 1, wx.ALL, 5)

        self.output_sizer = wx.BoxSizer(orient=wx.VERTICAL)
        self.output_box = wx.TextCtrl(self, style=wx.TE_MULTILINE)
        self.output_sizer.Add(self.output_box, 1, wx.ALL | wx.EXPAND, 5)

        self.panel_sizer.Add(self.file_picker_sizer, 0, wx.EXPAND, 5)
        self.panel_sizer.Add(self.run_sizer, 0, wx.EXPAND, 5)
        self.panel_sizer.Add(self.label_sizer, 0, wx.EXPAND, 5)
        self.panel_sizer.Add(self.output_sizer, 1, wx.EXPAND, 5)

        self.label_for_fasta.Show(False)
        self.label_box.Show(False)
        self.SetSizer(self.panel_sizer)
        self.Layout()

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