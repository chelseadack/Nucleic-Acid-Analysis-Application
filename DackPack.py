"""Nucleic Acid Analysis Application

Enter DNA or RNA strands, their concentrations, the temperature and salt conditions, then:
  1. Run Analysis   -> NUPACK calculates the equilibrium concentration of every complex (bar chart)
  2. Plot Structure -> NucDraw draws the MFE structure of a chosen complex, with base pairs
                       coloured by their equilibrium pair probability

Run with:  python3 GUI_run_script.py
"""
from __future__ import annotations

import sys

# import PySide6 before matplotlib, so matplotlib's Qt backend uses PySide6 even if PyQt is also installed
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QApplication, QComboBox, QDoubleSpinBox, QHBoxLayout, QLabel,
    QLineEdit, QMessageBox, QPushButton, QRadioButton, QSpinBox, QVBoxLayout, QWidget)

import matplotlib
matplotlib.use("QtAgg")

import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from nucdraw import NucDraw  # 10.5281/zenodo.15352138
from nupack import Model, SetSpec, Strand, Tube, tube_analysis

from GUI_logic import (M_TO_UM, MIN_CONCENTRATION_M, UM_TO_M, StrandInput, complex_label, complex_sort_key,
    convert_sequence, default_strand_name, dot_bracket_pairs, find_input_problems, has_only_valid_bases)

# ------------------------------
# PLOT SETTINGS
# ------------------------------
BASE_COLORS = {"DNA": dict(zip("ATCG", ["green", "red", "blue", "black"])),
               "RNA": dict(zip("AUCG", ["green", "red", "blue", "black"]))}
STRAND_COLORS = ["black", "blue", "red", "green", "orange", "darkviolet", "brown"]  # backbone colour of each strand
BASE_SIZES = {"big": (80, 6), "small": (30, 5)}  # marker area, letter font size
PAIR_CMAP = cm.jet

# ------------------------------
# STRAND INPUT ROW
# ------------------------------
class StrandRow(QWidget):
    """One row of the strand table: name, sequence, concentration and a remove button."""

    def __init__(self, parent_app: NupackApp):
        super().__init__()
        self.parent_app = parent_app

        layout = QHBoxLayout()
        layout.setSpacing(2)
        layout.setContentsMargins(0, 0, 0, 0)

        self.name = QLineEdit()
        self.name.setPlaceholderText("Name")
        self.name.setFixedWidth(70)

        self.sequence = QLineEdit()
        self.sequence.setPlaceholderText("e.g. ATGCGT")
        self.sequence.setMinimumWidth(140)

        self.conc = QDoubleSpinBox()
        self.conc.setValue(1.00)
        self.conc.setSuffix(" µM")
        self.conc.setFixedWidth(80)

        self.remove_btn = QPushButton("Remove")
        self.remove_btn.setProperty("danger", True)
        self.remove_btn.clicked.connect(self.remove_self)

        layout.addWidget(self.name)
        layout.addWidget(self.sequence)
        layout.addWidget(self.conc)
        layout.addWidget(self.remove_btn)

        self.setLayout(layout)

        self.sequence.textChanged.connect(self.validate_sequence)

    def validate_sequence(self) -> None:
        """Upper-case the sequence and outline it in red if it contains invalid bases."""
        cursor = self.sequence.cursorPosition()  # keep the cursor where the user is typing
        seq = self.sequence.text().upper()
        self.sequence.blockSignals(True)
        self.sequence.setText(seq)
        self.sequence.setCursorPosition(cursor)
        self.sequence.blockSignals(False)

        if has_only_valid_bases(seq, self.parent_app.material()):
            self.sequence.setStyleSheet("")
        else:
            self.sequence.setStyleSheet("border: 2px solid #e84118;")

    def convert_sequence(self, to_material: str) -> None:
        """Convert the sequence between DNA and RNA (T <-> U)."""
        self.sequence.setText(convert_sequence(self.sequence.text(), to_material))
        self.validate_sequence()

    def remove_self(self) -> None:
        self.parent_app.strand_rows.remove(self)
        self.deleteLater()

# ------------------------------
# MAIN APP
# ------------------------------
class NupackApp(QWidget):
    """Main window: model parameters, strand inputs and the structure viewer."""

    def __init__(self):
        super().__init__()

        self.setWindowTitle("Chelsea's NUPACK Desktop Tool")
        self.resize(480, 620)

        # results of the last analysis, used by Plot Structure
        self.tube_result_global = None
        self.complex_map = {}  # clean label -> NUPACK Complex
        self.run_material = None  # material & temperature of the last analysis, so the
        self.run_temperature = None  # structure title matches the analysis even if the inputs change

        self.eq_fig, self.eq_ax = None, None  # equilibrium concentration plot, reused between runs

        self.setStyleSheet(STYLESHEET)

        main_layout = QVBoxLayout()
        main_layout.setSpacing(2)
        main_layout.setContentsMargins(12, 12, 12, 12)
        self.setLayout(main_layout)

        main_layout.addWidget(self.section_label("1. Model Parameters"))
        main_layout.addWidget(self.model_card())

        main_layout.addWidget(self.section_label("2. Strand Sequences (5'-)"))
        main_layout.addWidget(self.strand_card())

        main_layout.addWidget(self.section_label("3. Structure Viewer"))
        main_layout.addWidget(self.structure_card())

    # ------------------------------
    # UI BUILDERS
    # ------------------------------
    def section_label(self, text: str) -> QLabel:
        label = QLabel(text)
        label.setProperty("section", True)
        return label

    def create_card(self, layout_inside: QVBoxLayout) -> QWidget:
        card = QWidget()
        card.setLayout(layout_inside)
        card.setObjectName("card")
        return card

    def form_row(self, label_text: str, widget: QWidget) -> QHBoxLayout:
        """Arrange a label and an input widget into a row."""
        row = QHBoxLayout()
        row.addWidget(QLabel(label_text))
        row.addStretch()
        row.addWidget(widget)
        return row

    def model_card(self) -> QWidget:
        layout = QVBoxLayout()
        layout.setSpacing(2)

        # Material
        row = QHBoxLayout()
        self.dna_radio = QRadioButton("DNA")
        self.rna_radio = QRadioButton("RNA")
        self.dna_radio.setChecked(True)
        row.addWidget(QLabel("Material:"))
        row.addWidget(self.dna_radio)
        row.addWidget(self.rna_radio)
        layout.addLayout(row)

        self.dna_radio.toggled.connect(self.on_material_change)

        # Temperature
        self.temperature = QDoubleSpinBox()
        self.temperature.setValue(37)
        layout.addLayout(self.form_row("Temperature (ºC):", self.temperature))

        # Sodium
        self.sodium = QDoubleSpinBox()
        self.sodium.setValue(1.0)
        self.sodium.setRange(0.05, 1.10)
        layout.addLayout(self.form_row("Sum of monovalent Na⁺ / K⁺ / NH₄⁺ ions (0.05 - 1.10 M):", self.sodium))

        # Magnesium
        self.magnesium = QDoubleSpinBox()
        self.magnesium.setValue(0.0)
        self.magnesium.setRange(0, 0.2)
        layout.addLayout(self.form_row("Sum of divalent Mg²⁺ ions (0 - 0.2 M):", self.magnesium))

        return self.create_card(layout)

    def strand_card(self) -> QWidget:
        layout = QVBoxLayout()

        self.strand_container = QVBoxLayout()
        self.strand_container.setSpacing(2)
        self.strand_rows = []

        self.add_strand()
        layout.addLayout(self.strand_container)

        self.add_strand_btn = QPushButton("Add Strand")
        self.add_strand_btn.setProperty("secondary", True)
        self.add_strand_btn.clicked.connect(self.add_strand)
        layout.addWidget(self.add_strand_btn)

        # Max size
        self.max_size = QSpinBox()
        self.max_size.setValue(4)
        layout.addLayout(self.form_row("Max complex size (strands):", self.max_size))

        # Run button
        self.run_btn = QPushButton("Run Analysis")
        self.run_btn.clicked.connect(self.run_analysis)
        layout.addWidget(self.run_btn)

        return self.create_card(layout)

    def structure_card(self) -> QWidget:
        layout = QVBoxLayout()
        layout.setSpacing(2)

        self.complex_dropdown = QComboBox()
        self.complex_dropdown.setFixedWidth(150)

        self.base_size = QComboBox()
        self.base_size.addItems(list(BASE_SIZES))
        self.base_size.setFixedWidth(80)

        self.rotate = QSpinBox()
        self.rotate.setRange(0, 180)
        self.rotate.setValue(80)
        self.rotate.setSuffix("º")

        layout.addLayout(self.form_row("Complex:", self.complex_dropdown))
        layout.addLayout(self.form_row("Base size:", self.base_size))
        layout.addLayout(self.form_row("Rotation:", self.rotate))

        self.plot_btn = QPushButton("Plot Structure")
        self.plot_btn.clicked.connect(self.plot_structure)
        layout.addWidget(self.plot_btn)

        return self.create_card(layout)

    # ------------------------------
    # INPUT
    # ------------------------------
    def material(self) -> str:
        return "RNA" if self.rna_radio.isChecked() else "DNA"

    def add_strand(self) -> None:
        row = StrandRow(self)
        self.strand_rows.append(row)
        self.strand_container.addWidget(row)

    def on_material_change(self) -> None:
        """Convert every sequence to the new material. RNA uses fixed salt conditions."""
        new_material = self.material()

        for row in self.strand_rows:
            row.convert_sequence(new_material)

        if new_material == "RNA":
            self.sodium.setValue(1.0)
            self.magnesium.setValue(0.0)
            self.sodium.setDisabled(True)
            self.magnesium.setDisabled(True)
        else:
            self.sodium.setDisabled(False)
            self.magnesium.setDisabled(False)

    def read_strands(self) -> list[StrandInput]:
        """Read the strands entered by the user, naming any unnamed strands."""
        return [StrandInput(name=row.name.text() or default_strand_name(i),
                            sequence=row.sequence.text(),
                            conc_um=row.conc.value())
                for i, row in enumerate(self.strand_rows)]

    # ------------------------------
    # ANALYSIS
    # ------------------------------
    def run_analysis(self) -> None:
        """Calculate the equilibrium concentration of every complex and plot them as a bar chart."""
        strands = self.read_strands()
        problems = find_input_problems(strands, self.material())
        if problems:
            QMessageBox.warning(self, "Invalid strand input", "\n".join(problems))
            return

        tube_strands = {Strand(s.sequence, name=s.name): s.conc_um * UM_TO_M for s in strands}
        t1 = Tube(strands=tube_strands, complexes=SetSpec(max_size=self.max_size.value()), name="Tube t1")
        model = Model(material=self.material(), ensemble="stacking", celsius=float(self.temperature.value()),
                      sodium=float(self.sodium.value()), magnesium=float(self.magnesium.value()))

        QApplication.setOverrideCursor(Qt.WaitCursor)  # NUPACK can take a few seconds for large simulations
        try:
            result = tube_analysis(tubes=[t1], compute=["pairs", "mfe"], model=model)
        except Exception as error:
            QMessageBox.critical(self, "NUPACK error", f"The analysis failed:\n{error}")
            return
        finally:
            QApplication.restoreOverrideCursor()

        # keep the complexes above the minimum concentration, sorted by size then name
        concs = result[t1].complex_concentrations
        filtered = {complex_label(k): (k, v) for k, v in concs.items() if v > MIN_CONCENTRATION_M}
        labels = sorted(filtered, key=complex_sort_key)
        values = [float(filtered[label][1]) * M_TO_UM for label in labels]

        if not labels:
            QMessageBox.information(self, "No complexes",
                f"No complex is above {MIN_CONCENTRATION_M * M_TO_UM:g} µM at equilibrium.")
            return

        self.plot_equilibrium(labels, values)

        # store the results for Plot Structure
        self.tube_result_global = result
        self.complex_map = {label: filtered[label][0] for label in labels}
        self.run_material = self.material()
        self.run_temperature = self.temperature.text()

        self.complex_dropdown.clear()
        self.complex_dropdown.addItems(labels)

    def plot_equilibrium(self, labels: list[str], values: list[float]) -> None:
        """Horizontal bar chart of equilibrium concentrations. The first run creates the figure,
        later runs redraw in the same window (or a new one, if it was closed)."""
        height = max(1.5, len(labels) * 0.5)
        if self.eq_fig is None or not plt.fignum_exists(self.eq_fig.number):
            self.eq_fig, self.eq_ax = plt.subplots(figsize=(10, height))
        else:
            self.eq_ax.clear()
            self.eq_fig.set_size_inches(10, height)  # fit the new number of complexes

        bars = self.eq_ax.barh(labels, values, color="red")
        self.eq_ax.invert_yaxis()
        self.eq_ax.set_xlabel("Equilibrium concentration (µM)")
        self.eq_ax.set_ylabel("Complex")
        self.eq_ax.set_xlim(0, max(values) + 0.5)

        # label the end of each bar with its concentration
        for bar, val in zip(bars, values):
            self.eq_ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height() / 2, f"{val:.2g} µM",
                            va="center", ha="left", fontsize=10)

        self.eq_fig.tight_layout()
        self.eq_fig.canvas.draw_idle()  # redraw the existing figure in place
        self.eq_fig.show()

    # ------------------------------
    # STRUCTURE PLOT
    # ------------------------------
    def plot_structure(self) -> None:
        """Draw the MFE structure of the chosen complex, with bases coloured by identity
        and base pairs coloured by their equilibrium probability."""
        name = self.complex_dropdown.currentText()
        if self.tube_result_global is None or name not in self.complex_map:
            QMessageBox.information(self, "No complex selected", "Run an analysis first, then choose a complex.")
            return

        complex_obj = self.complex_map[name]
        walker = self.tube_result_global[complex_obj]
        mfe = min(walker.mfe, key=lambda x: x.energy)
        structure = str(mfe.structure)

        # take the strands from the analysed complex (not the current inputs, which may have been edited since)
        strand_names = [strand.name for strand in complex_obj.strands]
        full_seq = "".join(str(strand) for strand in complex_obj.strands)
        if self.run_material == "RNA":
            full_seq = full_seq.replace("T", "U")

        # generate coordinates & draw the backbone and pairs
        nc = NucDraw(structure)
        nc.generate(degree=self.rotate.value())
        nc.plotter(12, bckwargs={"lw": 2, "color": "k"}, bpkwargs={"lw": 3, "c": "k"}, scwargs={"s": 0.01, "c": "k"})
        plt.title(f"MFE proxy structure for {self.run_material} complex {name} at {self.run_temperature}ºC "
                  f"is {mfe.energy:.2f} kcal/mol:")

        self.draw_bases(nc.coords, full_seq)
        nc.numbering_outside(20, 4, 10, {"fontsize": 10, "color": "k"})

        # if >1 strand, colour each backbone & add a legend
        if len(strand_names) > 1:
            strand_colors = STRAND_COLORS[:len(strand_names)]
            nc.multistrand_coloring(clr=strand_colors, bckwargs={"lw": 2})
            legend_elements = [Line2D([0], [0], color=color, lw=3, label=strand_name)
                               for color, strand_name in zip(strand_colors, strand_names)]
            plt.legend(handles=legend_elements, title="Strand:", loc="upper left")

        self.draw_pair_probabilities(nc.coords, structure, walker.pairs.to_array())

        plt.gcf().set_size_inches(7, 4)
        plt.show()

    def draw_bases(self, coords, full_seq: str) -> None:
        """Draw each base as a coloured circle with its letter."""
        base_colors = BASE_COLORS[self.run_material]
        base_diam, font_size = BASE_SIZES[self.base_size.currentText()]
        for (x, y), base in zip(coords, full_seq):
            plt.scatter(x, y, s=base_diam, facecolor=base_colors[base], edgecolor="none", zorder=3)
            plt.text(x - 0.3, y - 0.1, base, ha="center", va="center", fontsize=font_size, color="white", zorder=4)

    def draw_pair_probabilities(self, coords, structure: str, pair_probs) -> None:
        """Colour each MFE base pair by its equilibrium probability, with a colour bar."""
        ax = plt.gca()
        norm = mcolors.Normalize(vmin=0, vmax=1)
        for i, j in dot_bracket_pairs(structure).items():
            if i < j:
                (x1, y1), (x2, y2) = coords[i], coords[j]
                plt.plot([x1, x2], [y1, y2], color=PAIR_CMAP(norm(pair_probs[i, j])), linewidth=3, zorder=1)
        sm = cm.ScalarMappable(norm=norm, cmap=PAIR_CMAP)
        plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04, label="Probability of MFE state")

# ------------------------------
# STYLE
# ------------------------------
STYLESHEET = """
QWidget#card {background-color: white; border-radius: 10px; padding: 10px;}

QLabel[section="true"] {font-size: 16px; font-weight: bold; color: #2f3640;}

QLineEdit, QComboBox {background-color: white; border: 1px solid #dcdde1; border-radius: 6px; padding: 6px;}

QPushButton {background-color: #4078f2; color: white; border-radius: 6px; padding: 8px; font-weight: bold; border: none; background-clip: padding;}
QPushButton:hover {background-color: #2f5fd0;}
QPushButton:pressed {background-color: #244bb5;}

QPushButton[secondary="true"] {background-color: #dcdde1; color: #2f3640;}
QPushButton[secondary="true"]:hover {background-color: #c8cdd8;}
QPushButton[secondary="true"]:pressed {background-color: #b0b5c0;}

QPushButton[danger="true"] {background-color: #e84118;}
QPushButton[danger="true"]:hover {background-color: #c23616;}
QPushButton[danger="true"]:pressed {background-color: #a12a13;}

QSpinBox, QDoubleSpinBox {background-color: white; border: 1px solid #dcdde1; border-radius: 6px; padding: 4px; min-height: 14px;}
QSpinBox::up-button, QSpinBox::down-button, QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {width: 0px;}

QComboBox {border: 1px solid #dcdde1; border-radius: 6px; padding: 4px 6px; color: #2f3640;}
QComboBox QAbstractItemView {selection-background-color: #4078f2; selection-color: white; color: #2f3640;}
"""

# ------------------------------
# RUN
# ------------------------------
def main() -> None:
    app = QApplication(sys.argv)
    window = NupackApp()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
