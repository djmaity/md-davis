import collections
import os
import wx
from wx.lib.scrolledpanel import ScrolledPanel
from wx.lib.agw.floatspin import FloatSpin
from wx.lib.intctrl import IntCtrl
from wx.lib.filebrowsebutton import FileBrowseButton
import md_davis


class ElectrostaticsPanel(wx.lib.scrolledpanel.ScrolledPanel):

    def __init__(self, parent):
        ScrolledPanel.__init__(self, parent)

        self.panel_sizer = wx.BoxSizer(orient=wx.VERTICAL)

        self.name_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.panel_sizer.Add(self.name_sizer, 0, wx.EXPAND, 5)
        self.text_name = wx.StaticText(self, label='Name')
        self.name_sizer.Add(window=self.text_name, proportion=0, flag=wx.ALL,
                            border=5)
        self.name_textbox = wx.TextCtrl(self)
        self.name_sizer.Add(self.name_textbox, 1, wx.ALL, 5)

        self.output_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.panel_sizer.Add(self.output_sizer, 0, wx.EXPAND, 5)
        self.output_label = wx.StaticText(self, label='Output Directory')
        self.output_sizer.Add(window=self.output_label, proportion=0,
                              flag=wx.ALL, border=5)
        self.output_picker = wx.DirPickerCtrl(
            self,
            message='Select directory containing MD DaVis electrostatic files')
        self.output_sizer.Add(window=self.output_picker, proportion=1,
                              flag=wx.ALL, border=5)

        self.parameters_sizer = wx.BoxSizer(orient=wx.VERTICAL)
        self.panel_sizer.Add(self.parameters_sizer, 1, wx.EXPAND, 5)
        self.parameters_box = wx.TextCtrl(self, style=wx.TE_MULTILINE)
        self.parameters_sizer.Add(self.parameters_box, 1, wx.ALL | wx.EXPAND, 5)

        self.button_sizer = wx.BoxSizer(orient=wx.VERTICAL)
        self.panel_sizer.Add(self.button_sizer, 0, wx.EXPAND | wx.ALL, 5)
        self.button = wx.Button(self, label='Run DelPhi')
        self.button.Bind(wx.EVT_BUTTON, self.on_button_press)
        self.button_sizer.Add(self.button, proportion=0,
                              flag=wx.ALL | wx.CENTER, border=5)

        self.SetSizer(self.panel_sizer)
        self.Layout()
        self.Fit()

    def on_button_press(self, event):
        """
        Animate the electric field dynamics in PyMOL
        @param event: The event object
        """

        if self.output_picker.GetPath() or self.pymol_checkbox.GetValue():
            pass
        else:
            wx.MessageBox(
                "Please provide output file to save or select 'Show PyMOL Window'",
                "Message",
                wx.OK | wx.ICON_INFORMATION)
