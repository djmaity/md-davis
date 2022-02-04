import wx


class MyPanel(wx.Panel):

    def __init__(self, parent):
        super().__init__(parent)
        button = wx.Button(self, label='Get Sequence')
        button.Bind(wx.EVT_BUTTON, self.on_button_press)
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(button, proportion=0, flag=wx.ALL | wx.CENTER, border=5)
        self.SetSizer(main_sizer)

    def on_button_press(self, event):
        """
        Browse for a PDB file
        @param event: The event object
        """
        wildcard = "PDB files (*.pdb)|*.pdb"
        with wx.FileDialog(None, "Choose a file",
                           wildcard=wildcard,
                           style=wx.FD_OPEN) as dialog:
            if dialog.ShowModal() == wx.ID_OK:
                self.pdb_file.SetValue(dialog.GetPath())


class MyFrame(wx.Frame):

    def __init__(self):
        super().__init__(None, title='MD DaVis')
        panel = MyPanel(self)
        self.SetIcon(wx.Icon("MD_DaVis.ico"))
        self.Center()
        self.Show()


if __name__ == '__main__':
    app = wx.App()
    frame = MyFrame()
    app.MainLoop()
