import collections
import os
import wx
from wx.lib.scrolledpanel import ScrolledPanel
from wx.lib.intctrl import IntCtrl
from wx.lib.filebrowsebutton import FileBrowseButton
import md_davis


class MultipleFileBrowseButton(FileBrowseButton):

    def OnBrowse(self, event=None):
        """ Going to browse for file... """
        current = self.GetValue()
        directory = os.path.split(current)
        if os.path.isdir( current):
            directory = current
            current = ''
        elif directory and os.path.isdir( directory[0] ):
            current = directory[1]
            directory = directory [0]
        else:
            directory = self.startDirectory
            current = ''
        dlg = wx.FileDialog(self, self.dialogTitle, directory, current,
                            self.fileMask, self.fileMode)

        if dlg.ShowModal() == wx.ID_OK:
            self.SetValue(dlg.GetPaths())
        dlg.Destroy()

    def SetValue (self, value, callBack=1):
        """set current value of text control"""
        save = self.callCallback
        self.callCallback = callBack
        self.textControl.SetValue(' '.join([ '"'+_+'"' for _ in value]))
        self.callCallback =  save


class CollatePanel(wx.lib.scrolledpanel.ScrolledPanel):

    def __init__(self, parent):
        ScrolledPanel.__init__(self, parent)

        self.panel_sizer = wx.BoxSizer(orient=wx.VERTICAL)

        self.name_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.text_name = wx.StaticText(self, label="Name")
        self.name_sizer.Add(window=self.text_name, proportion=0, flag=wx.ALL, border=10)
        self.name_textbox = wx.TextCtrl(self)
        self.name_sizer.Add(self.name_textbox, 1, wx.ALL, 5)
        self.panel_sizer.Add(self.name_sizer, 0, wx.EXPAND, 5)

        self.label_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.text_label = wx.StaticText(self, label="Label")
        self.label_sizer.Add(window=self.text_label, proportion=0, flag=wx.ALL, border=10)
        self.label_textbox = wx.TextCtrl(self)
        self.label_sizer.Add(self.label_textbox, 1, wx.ALL, 5)
        self.panel_sizer.Add(self.label_sizer, 0, wx.EXPAND, 5)

        self.text_label_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.text_text_label = wx.StaticText(self, label="Text Label")
        self.text_label_sizer.Add(window=self.text_text_label, proportion=0, flag=wx.ALL, border=10)
        self.text_label_textbox = wx.TextCtrl(self)
        self.text_label_sizer.Add(self.text_label_textbox, 1, wx.ALL, 5)
        self.panel_sizer.Add(self.text_label_sizer, 0, wx.EXPAND, 5)

        self.structure_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.structure_label = wx.StaticText(self, label="Structure")
        self.structure_sizer.Add(window=self.structure_label, proportion=0, flag=wx.ALL, border=10)
        self.structure_picker = wx.FilePickerCtrl(self, path='', message="Select a structure",
                                                  wildcard='Structure files (*.pdb;*.gro)|*.pdb;*.gro',
                                                  style=wx.FLP_CHANGE_DIR | wx.FLP_USE_TEXTCTRL)
        self.structure_sizer.Add(window=self.structure_picker, proportion=1, flag=wx.ALL, border=5)
        self.panel_sizer.Add(self.structure_sizer, 0, wx.EXPAND, 5)

        self.trajectory_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.trajectory_label = wx.StaticText(self, label="Trajectory")
        self.trajectory_sizer.Add(window=self.trajectory_label, proportion=0, flag=wx.ALL, border=10)
        self.trajectory_picker = wx.FilePickerCtrl(self, path='', message="Select a structure",
                                                   wildcard='Structure files (*.trr;*.xtc)|*.trr;*.xtc',
                                                   style=wx.FLP_CHANGE_DIR | wx.FLP_USE_TEXTCTRL)
        self.trajectory_sizer.Add(window=self.trajectory_picker, proportion=1, flag=wx.ALL, border=5)
        self.panel_sizer.Add(self.trajectory_sizer, 0, wx.EXPAND, 5)

        # Time series: RMSD and Rg
        time_series_box = wx.StaticBox(self, label="Time series")
        self.time_series_sizer = wx.StaticBoxSizer(time_series_box, orient=wx.VERTICAL)
        self.panel_sizer.Add(self.time_series_sizer, 0, wx.EXPAND | wx.ALL, 5)

        self.rmsd_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.rmsd_label = wx.StaticText(self, label="RMSD")
        self.rmsd_sizer.Add(window=self.rmsd_label, proportion=0, flag=wx.ALL, border=10)
        self.rmsd_picker = wx.FilePickerCtrl(self, path='', message="Select a RMSD file",
                                             style=wx.FLP_CHANGE_DIR | wx.FLP_USE_TEXTCTRL)
        self.rmsd_sizer.Add(window=self.rmsd_picker, proportion=1, flag=wx.ALL, border=5)
        self.time_series_sizer.Add(self.rmsd_sizer, 0, wx.EXPAND, 5)

        self.rg_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.rg_label = wx.StaticText(self, label="Rg")
        self.rg_sizer.Add(window=self.rg_label, proportion=0, flag=wx.ALL, border=10)
        self.rg_picker = wx.FilePickerCtrl(self, path='', message="Select a structure",
                                           style=wx.FLP_CHANGE_DIR | wx.FLP_USE_TEXTCTRL)
        self.rg_sizer.Add(window=self.rg_picker, proportion=1, flag=wx.ALL, border=5)
        self.time_series_sizer.Add(self.rg_sizer, 0, wx.EXPAND, 5)

        # Dihedral Angles
        self.dihedral_box = wx.StaticBox(self)
        self.dihedral_sizer = wx.StaticBoxSizer(self.dihedral_box, orient=wx.VERTICAL)
        self.panel_sizer.Add(self.dihedral_sizer, 0, wx.EXPAND | wx.ALL, 5)

        self.calculate_dihedral = wx.CheckBox(self, label="Calculate dihedral angles")
        self.calculate_dihedral.Bind(wx.EVT_CHECKBOX, self.on_check_dihedral)
        self.dihedral_sizer.Add(window=self.calculate_dihedral, proportion=0, flag=wx.ALL, border=10)

        self.chunk_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.chunk_text = wx.StaticText(self, label="Chunk Size")
        self.chunk_sizer.Add(window=self.chunk_text, proportion=0, flag=wx.ALL, border=10)
        self.chunk_intbox = IntCtrl(self, value=1000)
        self.chunk_sizer.Add(self.chunk_intbox, 1, wx.ALL, 5)
        self.dihedral_sizer.Add(self.chunk_sizer, 0, wx.EXPAND, 5)
        self.chunk_text.Show(False)
        self.chunk_intbox.Show(False)

        # Residue Properties
        self.residue_box = wx.StaticBox(self, label='Residue properties')
        self.residue_sizer = wx.StaticBoxSizer(self.residue_box, orient=wx.VERTICAL)
        self.panel_sizer.Add(self.residue_sizer, 0, wx.EXPAND | wx.ALL, 5)

        self.secstr_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.secstr_label = wx.StaticText(self, label="DSSP")
        self.secstr_sizer.Add(window=self.secstr_label, proportion=0, flag=wx.ALL, border=10)
        self.secstr_picker = wx.FilePickerCtrl(self, path='', message="DSSP .dat file",
                                               wildcard='DSSP files (*.dat)|*.dat',
                                               style=wx.FLP_CHANGE_DIR | wx.FLP_USE_TEXTCTRL)
        self.secstr_sizer.Add(window=self.secstr_picker, proportion=1, flag=wx.ALL, border=5)
        self.residue_sizer.Add(self.secstr_sizer, 0, wx.EXPAND, 5)

        self.rmsf_box = wx.StaticBox(self, label='Root-mean-square fluctuation (RMSF)')
        self.rmsf_box_sizer = wx.StaticBoxSizer(self.rmsf_box, orient=wx.VERTICAL)
        self.residue_sizer.Add(self.rmsf_box_sizer, 0, wx.EXPAND | wx.ALL, 5)

        # TODO: Allow multiple file input to RMSF
        self.rmsf_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.rmsf_label = wx.StaticText(self, label="RMSF")
        self.rmsf_sizer.Add(window=self.rmsf_label, proportion=0, flag=wx.ALL, border=10)
        self.rmsf_picker = wx.FilePickerCtrl(self, path='', message="RMSF file",
                                             style=wx.FLP_CHANGE_DIR | wx.FLP_USE_TEXTCTRL)
        # self.rmsf_picker = MultipleFileBrowseButton(self,  dialogTitle="Select RMSF file", labelText='RMSF',
        #                                     fileMode=wx.FD_MULTIPLE)
        self.rmsf_sizer.Add(window=self.rmsf_picker, proportion=1, flag=wx.ALL, border=5)
        self.rmsf_box_sizer.Add(self.rmsf_sizer, 0, wx.EXPAND, 5)

        self.sasa_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.sasa_label = wx.StaticText(self, label="SASA")
        self.sasa_sizer.Add(window=self.sasa_label, proportion=0, flag=wx.ALL, border=10)
        self.sasa_picker = wx.FilePickerCtrl(self, path='', message="Select a SASA file",
                                             style=wx.FLP_CHANGE_DIR | wx.FLP_USE_TEXTCTRL)
        self.sasa_sizer.Add(window=self.sasa_picker, proportion=1, flag=wx.ALL, border=5)
        self.residue_sizer.Add(self.sasa_sizer, 0, wx.EXPAND, 5)

        self.electrostatics_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.electrostatics_label = wx.StaticText(self, label="Electrostatics")
        self.electrostatics_sizer.Add(window=self.electrostatics_label, proportion=0, flag=wx.ALL, border=10)
        self.electrostatics_picker = wx.DirPickerCtrl(
            self, message="Select directory containing MD DaVis electrostatics files")
        self.electrostatics_sizer.Add(window=self.electrostatics_picker, proportion=1, flag=wx.ALL, border=5)
        self.residue_sizer.Add(self.electrostatics_sizer, 0, wx.EXPAND, 5)

        self.output_sizer = wx.BoxSizer(orient=wx.HORIZONTAL)
        self.output_label = wx.StaticText(self, label="HDF File")
        self.output_sizer.Add(window=self.output_label, proportion=0, flag=wx.ALL, border=10)
        self.output_picker = wx.FilePickerCtrl(
            self, path='', message="Output HDF file",
            wildcard='HDF file (*.hdf;*.h5;*.hdf5;*.he5)|*.h5;*.hdf5;*.he5;*.hdf',
            style=wx.FLP_SAVE | wx.FLP_CHANGE_DIR | wx.FLP_USE_TEXTCTRL)
        self.output_sizer.Add(window=self.output_picker, proportion=1, flag=wx.ALL, border=5)
        self.panel_sizer.Add(self.output_sizer, 0, wx.EXPAND, 5)

        self.button_sizer = wx.BoxSizer(orient=wx.VERTICAL)
        self.button = wx.Button(self, label='Create HDF File')
        self.button.Bind(wx.EVT_BUTTON, self.on_button_press)
        self.button_sizer.Add(self.button, proportion=0, flag=wx.ALL | wx.CENTER, border=5)
        self.panel_sizer.Add(self.button_sizer, 0, wx.EXPAND | wx.ALL, 5)

        self.SetSizer(self.panel_sizer)
        self.Layout()
        self.SetupScrolling()

    def on_check_dihedral(self, event):
        """
        Add chunk box on selecting to calculate dihedral angles
        @param event: The event object
        """
        if self.calculate_dihedral.GetValue():
            self.chunk_text.Show()
            self.chunk_intbox.Show()
            self.Layout()
            self.SetupScrolling(scrollToTop=False)
        else:
            self.chunk_text.Show(False)
            self.chunk_intbox.Show(False)
            self.Layout()
            self.SetupScrolling(scrollToTop=False)


    def on_button_press(self, event):
        """
        Browse for a PDB file
        @param event: The event object
        """
        data = collections.defaultdict(dict)
        if self.name_textbox.GetValue():
            data['name'] = self.name_textbox.GetValue()
        if self.label_textbox.GetValue():
            data['label'] = self.label_textbox.GetValue()
        if self.text_label_textbox.GetValue():
            data['text_label'] = self.text_label_textbox.GetValue()
        if self.output_picker.GetPath():
            data['output'] = self.output_picker.GetPath()
        if self.structure_picker.GetPath():
            data['structure'] = self.structure_picker.GetPath()
        if self.trajectory_picker.GetPath():
            data['trajectory'] = self.trajectory_picker.GetPath()
        if self.rmsd_picker.GetPath():
            data['timeseries']['rmsd'] = self.rmsd_picker.GetPath()
        if self.rg_picker.GetPath():
            data['timeseries']['rg'] = self.rg_picker.GetPath()
        if self.secstr_picker.GetPath():
            data['residue_property']['secondary_structure'] = self.secstr_picker.GetPath()
        if self.rmsf_picker.GetPath():
            data['residue_property']['rmsf'] = self.rmsf_picker.GetPath()
        if self.sasa_picker.GetPath():
            data['residue_property']['sasa'] = self.sasa_picker.GetPath()
        if self.electrostatics_picker.GetPath():
            data['residue_property']['electrostatics'] = self.electrostatics_picker.GetPath()
        if  self.calculate_dihedral.GetValue():
            data['dihedral']['chunk'] = self.chunk_intbox.GetValue()

        if md_davis.collate.create_hdf(data):
            print('Collated all data into ' + self.output_picker.GetPath())

