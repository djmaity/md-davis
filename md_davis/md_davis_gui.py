import wx
from md_davis.gui.sequence_panel import SequencePanel
from md_davis.gui.collate_panel import CollatePanel
from md_davis.gui.landscape_panel import LandscapePanel
# from md_davis.gui.residue_panel import ResiduePanel


class MainPanel(wx.Panel):

    def __init__(self, parent):
        super().__init__(parent)

        main_sizer = wx.BoxSizer(orient=wx.VERTICAL)
        notebook = wx.Notebook(self)
        notebook.AddPage(LandscapePanel(notebook), 'Landscape')
        notebook.AddPage(CollatePanel(notebook), 'Collate')
        notebook.AddPage(SequencePanel(notebook), 'Sequence')
        # notebook.AddPage(ResiduePanel(notebook), 'Residue')
        main_sizer.Add(notebook, 1, wx.ALL | wx.EXPAND, 5)
        self.SetSizer(main_sizer)
        self.Layout()


class MainFrame(wx.Frame):

    def __init__(self):
        super().__init__(None, title='MD DaVis')
        self.SetIcon(wx.Icon('MD_DaVis.ico'))
        self.panel = MainPanel(self)

        # # Status Bar
        # self.CreateStatusBar()
        #
        # # File Menu
        # self.File = wx.Menu()
        # self.load = wx.MenuItem(self.File, wx.ID_ANY, 'Load', 'Load a input parameter file', wx.ITEM_NORMAL)
        # self.File.Append(self.load)
        #
        # # Tools Menu
        # self.Tools = wx.Menu()
        # self.sequence = wx.MenuItem(self.Tools, wx.ID_ANY, 'Get Sequence', 'Get the Sequence from a PDB file', wx.ITEM_NORMAL)
        # self.Tools.Append(self.sequence)
        # # self.Bind(wx.EVT_MENU, self.OnSequenceMenu, menuItem)
        #
        # self.collate = wx.MenuItem(self.Tools, wx.ID_ANY, 'Collate Data', wx.EmptyString, wx.ITEM_NORMAL)
        # self.Tools.Append(self.collate)
        #
        # self.landscape = wx.MenuItem(self.Tools, wx.ID_ANY, 'Free Energy Landscape', wx.EmptyString, wx.ITEM_NORMAL)
        # self.Tools.Append(self.landscape)
        #
        # self.residue = wx.MenuItem(self.Tools, wx.ID_ANY, 'Residue Properties', wx.EmptyString, wx.ITEM_NORMAL)
        # self.Tools.Append(self.residue)
        #
        # self.Help = wx.Menu()
        # self.about = wx.MenuItem(self.Help, wx.ID_ANY, 'About', wx.EmptyString, wx.ITEM_NORMAL)
        # self.Tools.Append(self.about)
        #
        # # Create and Add Menus to Menu bar
        # self.menubar = wx.MenuBar()
        # self.menubar.Append(self.Tools, 'Tools')
        # self.menubar.Append(self.Help, 'Help')
        #
        # self.SetMenuBar(self.menubar)

        self.SetMinSize(wx.Size(width=400, height=400))
        self.Center()
        self.Show()


if __name__ == '__main__':
    app = wx.App()
    frame = MainFrame()
    app.MainLoop()
