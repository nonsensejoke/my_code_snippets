'''
Optimize
Described at PyMOL wiki: http://www.pymolwiki.org/index.php/optimize

Author : Osvaldo Martin
email: aloctavodia@gmail.com
Date: august 2014
License: MIT License
Version 0.9

how to use: 
 - command line 
    minimize loop, forcefield=GAFF, nsteps0=100, cutoff=1,update_in_place=1 

'''

import sys
from pymol import cmd

try:
    from openbabel import openbabel as ob
except ImportError:
    print('<' * 80 + '''

Optimize plug-in needs openbabel to be installed in your system, please follow the instructions at
http://openbabel.org/wiki/Get_Open_Babel

''' + '>' * 80)

from pymol.Qt import QtCore, QtWidgets, QtGui
Qt = QtCore.Qt


def __init_plugin__(app=None) -> None:
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('OpenBabel Optimize MOD', run_plugin_gui)

def run_plugin_gui() -> None:
    """
    Creates and shows the Optimize dialog.
    """
    global dialog

    if dialog is None:
        dialog = make_dialog()

    dialog.show()

dialog = None


def make_dialog():
    """
    Create the dialog
    """
    return OptimizeDialog()

class OptimizeDialog(QtWidgets.QDialog):
    """
    Dialog for optimization settings.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(" Optimize ")
        self.setMinimumWidth(450)

        # Main vertical layout
        main_layout = QtWidgets.QVBoxLayout(self)

        # Header label (replaces the original Tkinter Label)
        header_label = QtWidgets.QLabel("Optimize: Let's find that minimum!")
        header_label.setStyleSheet("background-color: black; color: white; padding: 4px;")
        header_label.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(header_label)

        # Tab widget (replaces Pmw.NoteBook)
        self.tabs = QtWidgets.QTabWidget()
        main_layout.addWidget(self.tabs)

        # Create and add the tabs
        self._create_local_opt_tab()
        self._create_global_opt_tab()
        self._create_about_tab()

    def _create_local_opt_tab(self):
        """
        Creates the 'Local optimization' tab."""
        local_tab = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(local_tab)
        layout.setContentsMargins(10, 10, 10, 10)

        # Group box for options (replaces Pmw.Group)
        group_box = QtWidgets.QGroupBox("Minimization options")
        grid = QtWidgets.QGridLayout(group_box) # Using a grid for clean alignment

        # Get current molecule selections from PyMOL
        pymol_selections = cmd.get_names('all')
        if not pymol_selections:
            pymol_selections = ['all']

        grid.addWidget(QtWidgets.QLabel("Force Field:"), 0, 0)
        self.local_ff_combo = QtWidgets.QComboBox()
        self.local_ff_combo.addItems(['GAFF', 'MMFF94s', 'MMFF94', 'UFF', 'Ghemical'])
        self.local_ff_combo.setCurrentText('MMFF94s')
        grid.addWidget(self.local_ff_combo, 0, 1)

        # Method
        grid.addWidget(QtWidgets.QLabel("Method:"), 1, 0)
        self.local_method_combo = QtWidgets.QComboBox()
        self.local_method_combo.addItems(['Conjugate Gradients', 'Steepest Descent'])
        grid.addWidget(self.local_method_combo, 1, 1)

        # Steps (replaces Label and Entry with QLabel and QLineEdit)
        grid.addWidget(QtWidgets.QLabel("Steps:"), 2, 0)
        self.entry_nsteps0 = QtWidgets.QLineEdit("500")
        grid.addWidget(self.entry_nsteps0, 2, 1)

        # Convergence
        grid.addWidget(QtWidgets.QLabel("Convergence:"), 3, 0)
        self.entry_conv = QtWidgets.QLineEdit("0.0001")
        grid.addWidget(self.entry_conv, 3, 1)

        # Selection
        grid.addWidget(QtWidgets.QLabel("Selection:"), 4, 0)
        self.sel0_value = QtWidgets.QLineEdit(pymol_selections[0])
        grid.addWidget(self.sel0_value, 4, 1)

        # Cutoff Radio Buttons (replaces Tkinter.Radiobutton)
        self.no_cutoff_radio = QtWidgets.QRadioButton("No cutoff")
        self.no_cutoff_radio.setChecked(True)
        grid.addWidget(self.no_cutoff_radio, 5, 0, 1, 2)

        self.use_cutoff_radio = QtWidgets.QRadioButton("Use cutoff")
        grid.addWidget(self.use_cutoff_radio, 6, 0, 1, 2)

        # Van der Waals Cutoff
        grid.addWidget(QtWidgets.QLabel("Van der Waals Cutoff (Å):"), 7, 0)
        self.entry_vdw = QtWidgets.QLineEdit("6.0")
        self.entry_vdw.setEnabled(False)
        grid.addWidget(self.entry_vdw, 7, 1)

        # Electrostatic Cutoff
        grid.addWidget(QtWidgets.QLabel("Electrostatic Cutoff (Å):"), 8, 0)
        self.entry_elec = QtWidgets.QLineEdit("8.0")
        self.entry_elec.setEnabled(False)
        grid.addWidget(self.entry_elec, 8, 1)

        # Connect signal to slot for enabling/disabling cutoff entries
        self.use_cutoff_radio.toggled.connect(self._toggle_cutoff_entries)

        # Update in place option
        self.update_in_place_check = QtWidgets.QCheckBox(
            "Update coordinates in place (no new object)")
        self.update_in_place_check.setToolTip(
            "When checked, optimized coordinates are written back to the original\n"
            "selection via cmd.load_coords instead of creating a new PyMOL object.\n"
            "Preserves colors, representations, and selections on the parent object.")
        grid.addWidget(self.update_in_place_check, 9, 0, 1, 2)

        layout.addWidget(group_box)
        layout.addStretch()

        # Minimize Button (replaces Tkinter.Button)
        minimize_button = QtWidgets.QPushButton("Minimize")
        minimize_button.clicked.connect(self._set_minimize)
        layout.addWidget(minimize_button, alignment=Qt.AlignRight)

        self.tabs.addTab(local_tab, " Local optimization ")

    def _create_global_opt_tab(self):
        """
        Creates the 'Global Optimization' tab.
        """
        global_tab = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(global_tab)
        layout.setContentsMargins(10, 10, 10, 10)

        group_box = QtWidgets.QGroupBox("Conformational Search options")
        grid = QtWidgets.QGridLayout(group_box)

        pymol_selections = cmd.get_names('all')
        if not pymol_selections:
            pymol_selections = ['all']

        # --- Widgets for Global Optimization ---
        # Force Field
        grid.addWidget(QtWidgets.QLabel("Force Field:"), 0, 0)
        self.global_ff_combo = QtWidgets.QComboBox()
        self.global_ff_combo.addItems(['GAFF', 'MMFF94s', 'MMFF94', 'UFF', 'Ghemical'])
        self.global_ff_combo.setCurrentText('MMFF94s')
        grid.addWidget(self.global_ff_combo, 0, 1)

        # Method
        grid.addWidget(QtWidgets.QLabel("Method:"), 1, 0)
        self.conf_method_combo = QtWidgets.QComboBox()
        self.conf_method_combo.addItems(['Weighted', 'Random', 'Systematic'])
        self.conf_method_combo.currentTextChanged.connect(self._toggle_conformer_entries)
        grid.addWidget(self.conf_method_combo, 1, 1)

        # Steps
        grid.addWidget(QtWidgets.QLabel("Steps (for optimization):"), 2, 0)
        self.entry_nsteps1 = QtWidgets.QLineEdit("500")
        grid.addWidget(self.entry_nsteps1, 2, 1)

        # Conformers
        grid.addWidget(QtWidgets.QLabel("Conformers (to generate):"), 3, 0)
        self.entry_conformers = QtWidgets.QLineEdit("25")
        grid.addWidget(self.entry_conformers, 3, 1)

        # Lowest Energy Conformers
        grid.addWidget(QtWidgets.QLabel("Lowest Energy (to keep):"), 4, 0)
        self.entry_lowest = QtWidgets.QLineEdit("5")
        grid.addWidget(self.entry_lowest, 4, 1)

        # Selection
        grid.addWidget(QtWidgets.QLabel("Selection:"), 5, 0)
        self.sel1_value = QtWidgets.QLineEdit(pymol_selections[0])
        grid.addWidget(self.sel1_value, 5, 1)

        layout.addWidget(group_box)
        layout.addStretch()

        # Search Button
        search_button = QtWidgets.QPushButton("Search")
        search_button.clicked.connect(self._set_conf_search)
        layout.addWidget(search_button, alignment=Qt.AlignRight)

        self.tabs.addTab(global_tab, " Global Optimization ")
        self._toggle_conformer_entries(self.conf_method_combo.currentText())

    def _create_about_tab(self):
        """
        Creates the 'About' tab.
        """
        about_tab = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(about_tab)
        layout.setContentsMargins(15, 15, 15, 15)

        about_text = """
Optimize provides a PyMOL graphical interface to some
of the many options available in openbabel (openbabel.org).


If you find Optimize useful great! 
If you don't and have some suggestions or comments 
to do please write to me (aloctavodia@gmail.com).
"""
        label = QtWidgets.QLabel(about_text)
        label.setAlignment(Qt.AlignCenter)
        label.setWordWrap(True)

        layout.addWidget(label)
        layout.addStretch()

        self.tabs.addTab(about_tab, "    About   ")

    def _toggle_cutoff_entries(self, checked):
        """
        Enables or disables cutoff fields based on the radio button state.
        """
        self.entry_vdw.setEnabled(checked)
        self.entry_elec.setEnabled(checked)

    def _toggle_conformer_entries(self, method_text):
        """
        Enables or disables conformer fields based on the selected method.
        """
        is_systematic = (method_text == 'Systematic')
        self.entry_conformers.setEnabled(not is_systematic)
        self.entry_lowest.setEnabled(not is_systematic)

    def _set_minimize(self):
        """
        Reads values from the GUI and calls the minimize() function.
        """
        try:
            forcefield = self.local_ff_combo.currentText()
            method = self.local_method_combo.currentText()
            nsteps0 = int(self.entry_nsteps0.text())
            conv = float(self.entry_conv.text())
            cutoff = self.use_cutoff_radio.isChecked()
            cut_vdw = float(self.entry_vdw.text())
            cut_elec = float(self.entry_elec.text())
            selection = self.sel0_value.text()
            update_in_place = self.update_in_place_check.isChecked()

            print("Starting minimization...")
            minimize(selection, forcefield, method, nsteps0, conv, cutoff,
                     cut_vdw, cut_elec, update_in_place=update_in_place)
        except ValueError as e:
            print(f"Error: Invalid input value. Please check numbers. Details: {e}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

    def _set_conf_search(self):
        """
        Reads values from the GUI and calls the conf_search() function.
        """
        try:
            forcefield = self.global_ff_combo.currentText()
            conf_method = self.conf_method_combo.currentText()
            nsteps1 = int(self.entry_nsteps1.text())
            conformers = int(self.entry_conformers.text())
            lowest_conf = int(self.entry_lowest.text())
            selection = self.sel1_value.text()

            print("Starting conformational search...")
            conf_search(selection, forcefield, conf_method, nsteps1, conformers, lowest_conf)
        except ValueError as e:
            print(f"Error: Invalid input value. Please check numbers. Details: {e}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")


def minimize(selection: str ='all', forcefield: str ='MMFF94s',
             method: str ='Conjugate Gradients', nsteps0: int = 500,
             conv: float = 0.0001, cutoff: bool = False,
             cut_vdw: float = 6.0, cut_elec: float = 8.0,
             update_in_place: bool = False) -> None:
    """
DESCRIPTION

    Minimize the energy of a molecule using OpenBabel's force fields.

ARGUMENTS

    selection = string: The selection string for the molecule to minimize.

    forcefield = string: The force field to use (e.g., 'GAFF', 'MMFF94s', 'UFF', 'Ghemical').

    method = string: The optimization method ('Conjugate Gradients' or 'Steepest Descent').

    nsteps0 = int: The number of optimization steps.

    conv = float: The convergence criterion.

    cutoff = bool: Whether to use cutoff for van der Waals and electrostatic interactions.

    cut_vdw = float: The cutoff distance for van der Waals interactions.

    cut_elec = float: The cutoff distance for electrostatic interactions.

    update_in_place = bool: When True, write optimized coordinates back into the original
        selection via cmd.load_coords instead of deleting and reloading as a new object.
        Preserves colors, representations, and selections on the parent object.

SEE ALSO

    conf_search

    """
    # cmd.extend passes all CLI arguments as strings — coerce to correct types upfront
    nsteps0 = int(nsteps0)
    conv = float(conv)
    cut_vdw = float(cut_vdw)
    cut_elec = float(cut_elec)
    cutoff = bool(int(cutoff)) if isinstance(cutoff, str) else bool(cutoff)
    update_in_place = bool(int(update_in_place)) if isinstance(update_in_place, str) else bool(update_in_place)

    print('[DEBUG] === minimize() called ===')
    print(f'[DEBUG] selection="{selection}", forcefield={forcefield}, method={method}')
    print(f'[DEBUG] nsteps0={nsteps0} (int), conv={conv} (float), cutoff={cutoff} (bool), update_in_place={update_in_place} (bool)')

    # --- Step 1: Validate selection and export from PyMOL ---
    pymol_atom_count = cmd.count_atoms(selection)
    print(f'[DEBUG] PyMOL atom count for "{selection}": {pymol_atom_count}')
    if pymol_atom_count == 0:
        print(f'[DEBUG] ERROR: selection "{selection}" matched 0 atoms. Aborting.')
        return

    mol_string = cmd.get_str('mol', selection)
    name = cmd.get_legal_name(selection)
    print(f'[DEBUG] Legal name resolved to: "{name}"')
    print(f'[DEBUG] Exported MOL string length: {len(mol_string)} chars')

    # --- Step 2: Parse with OpenBabel ---
    obconversion = ob.OBConversion()
    obconversion.SetInAndOutFormats('mol', 'mol')
    mol = ob.OBMol()
    read_ok = obconversion.ReadString(mol, mol_string)
    ob_atom_count = mol.NumAtoms()
    print(f'[DEBUG] OBConversion.ReadString success: {read_ok}')
    print(f'[DEBUG] OBMol atom count after parsing: {ob_atom_count}')

    if ob_atom_count == 0:
        print('[DEBUG] ERROR: OBMol has 0 atoms. The MOL string may be malformed. Aborting.')
        return

    if ob_atom_count != pymol_atom_count:
        print(f'[DEBUG] WARNING: atom count mismatch — PyMOL={pymol_atom_count}, OB={ob_atom_count}')
        print('[DEBUG] This may indicate H atoms were added/removed during MOL parsing.')

    # --- Step 3: Force field setup ---
    ff = ob.OBForceField.FindForceField(forcefield)
    if ff is None:
        print(f'[DEBUG] ERROR: Force field "{forcefield}" not found by OpenBabel. Aborting.')
        return
    ff_ok = ff.Setup(mol)
    print(f'[DEBUG] ForceField.Setup() returned: {ff_ok}')
    if not ff_ok:
        print(f'[DEBUG] WARNING: Force field setup failed (possibly missing atom type parameters).')
        print(f'[DEBUG] Results may be incorrect. Consider switching to MMFF94s or UFF.')

    # --- Step 4: Optional cutoff ---
    if cutoff:
        ff.EnableCutOff(True)
        ff.SetVDWCutOff(cut_vdw)
        ff.SetElectrostaticCutOff(cut_elec)
        print(f'[DEBUG] Cutoff enabled — VDW: {cut_vdw} Å, Electrostatic: {cut_elec} Å')

    # --- Step 5: Run optimization ---
    print(f'[DEBUG] Starting {method} ({nsteps0} steps, convergence={conv}) ...')
    if method == 'Conjugate Gradients':
        ff.ConjugateGradients(nsteps0, conv)
    else:
        ff.SteepestDescent(nsteps0, conv)

    ff.GetCoordinates(mol)
    nrg = ff.Energy()
    print(f'[DEBUG] Optimization complete. Final energy: {nrg:.4f} {ff.GetUnit()}')

    # --- Step 6: Write results back to PyMOL ---
    if update_in_place:
        # Safety check: atom counts must match for index-based coordinate mapping
        if pymol_atom_count != ob_atom_count:
            print(f'[DEBUG] ERROR: cannot update in place — atom count mismatch '
                  f'(PyMOL={pymol_atom_count}, OB={ob_atom_count}).')
            print('[DEBUG] OpenBabel may have added/removed atoms during setup. Aborting.')
            return

        coords = [[mol.GetAtom(i).GetX(), mol.GetAtom(i).GetY(), mol.GetAtom(i).GetZ()]
                  for i in range(1, ob_atom_count + 1)]
        print(f'[DEBUG] Extracted {len(coords)} coordinate sets from OBMol')
        print(f'[DEBUG] Sample — atom 1: {coords[0]}, atom -1: {coords[-1]}')
        print(f'[DEBUG] Calling cmd.load_coords on "{selection}", state=1 ...')
        cmd.load_coords(coords, selection, state=1)
        print(f'[DEBUG] cmd.load_coords completed. Original object preserved.')
        display_name = selection
    else:
        mol_string = obconversion.WriteString(mol)
        print(f'[DEBUG] Written MOL string length: {len(mol_string)} chars')
        cmd.delete(name)
        if name == 'all':
            name = 'all_'
        cmd.read_molstr(mol_string, name, state=0, finish=1, discrete=1)
        print(f'[DEBUG] Loaded new object "{name}" into PyMOL.')
        display_name = name

    print('#########################################')
    print('The Energy of %s is %8.2f %s       ' % (display_name, nrg, ff.GetUnit()))
    print('#########################################')


def conf_search(selection: str = 'all', forcefield: str = 'MMFF94s',
                method: str = 'Weighted', nsteps1: int = 500,
                conformers: int = 25, lowest_conf: int = 5):
    """
DESCRIPTION

    Perform a conformational search on a molecule using OpenBabel's force fields.

ARGUMENTS

    selection = string: The selection string for the molecule to search.

    forcefield = string: The force field to use (e.g., 'GAFF', 'MMFF94s', 'UFF', 'Ghemical').

    method = string: The search method for global-minimum ('Weighted', 'Random', or 'Systematic').

    nsteps1 = int: The number of optimization steps for each conformer.

    conformers = int: The number of conformers to be analyzed.

    lowest_conf = int: The number of lowest energy conformers to keep.

SEE ALSO

    minimize
    """
    mol_string = cmd.get_str('mol', selection)
    name = cmd.get_legal_name(selection)
    obconversion = ob.OBConversion()
    obconversion.SetInAndOutFormats('mol', 'mol')
    mol = ob.OBMol()
    obconversion.ReadString(mol, mol_string)
    ff = ob.OBForceField.FindForceField(forcefield) ## GAFF, MMFF94s, MMFF94, UFF, Ghemical
    ff.Setup(mol)
    if method == 'Weighted':
        ff.WeightedRotorSearch(conformers, nsteps1)
    elif method == 'Random':
        ff.RandomRotorSearch(conformers, nsteps1)
    else:
        ff.SystematicRotorSearch(nsteps1)
    if name == 'all':
        name = 'all_'
    if method in ['Weighted', 'Random']:
        ff.GetConformers(mol)
        print('##############################################')
        print('   Conformer    |         Energy      |  RMSD')
        nrg_unit = ff.GetUnit()
        rmsd = 0
        ff.GetCoordinates(mol)
        nrg = ff.Energy()
        conf_list = []
        for i in range(conformers):
            mol.SetConformer(i) 
            ff.Setup(mol)
            nrg = ff.Energy()
            conf_list.append((nrg, i))
        conf_list.sort()
        lenght_conf_list = len(conf_list)
        if lowest_conf > lenght_conf_list:
            lowest_conf = lenght_conf_list
        for i in range(lowest_conf):
            nrg, orden = conf_list[i]
            name_n = '%s%02d' % (name, i)
            cmd.delete(name_n)
            mol.SetConformer(orden) 
            mol_string = obconversion.WriteString(mol)
            cmd.read_molstr(mol_string, name_n,state=0,finish=1,discrete=1)
            if i != 0:
                rmsd = cmd.fit(name_n, '%s00' % name, quiet=1)
            print('%15s | %10.2f%9s |%6.1f'    % (name_n, nrg, nrg_unit, rmsd))
        print('##############################################')
    else:
        ff.GetCoordinates(mol)
        nrg = ff.Energy()
        mol_string = obconversion.WriteString(mol)
        cmd.delete(name)
        cmd.read_molstr(mol_string, name,state=0,finish=1,discrete=1)
        print('#########################################')
        print('The Energy of %s is %8.2f %s       '  % (name, nrg, ff.GetUnit()))
        print('#########################################')

cmd.extend('minimize', minimize)
cmd.extend('conf_search', conf_search)
