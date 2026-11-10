import asyncio
from dataclasses import asdict
import html
import json

import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
from IPython.display import display

from cidep import Config, Group, Radical, Shape, broaden, simulate
from storage import ARRAYS, import_legacy, save, validate


class RadicalEditor:
    def __init__(self, radical):
        self.g = widgets.FloatText(value=radical.g, description="g-value")
        self.rows = widgets.VBox()
        self.add = widgets.Button(icon="plus", tooltip="Add nuclear group", layout=widgets.Layout(width="36px"))
        self.add.on_click(lambda _: self.append(Group()))
        self.widget = widgets.VBox([self.g, widgets.HBox([
            widgets.Label("Nuclei", layout=widgets.Layout(width="75px")),
            widgets.Label("Coupling / G", layout=widgets.Layout(width="110px")),
            widgets.Label("Spin", layout=widgets.Layout(width="75px"))]), self.rows, self.add])
        self.set(radical)

    def append(self, group):
        count = widgets.BoundedIntText(value=group.count, min=1, max=15, layout=widgets.Layout(width="75px"))
        coupling = widgets.FloatText(value=group.coupling, layout=widgets.Layout(width="110px"))
        spin = widgets.BoundedFloatText(value=group.spin, min=0, max=5, step=0.5, layout=widgets.Layout(width="75px"))
        remove = widgets.Button(icon="trash", tooltip="Remove group", layout=widgets.Layout(width="36px"))
        row = widgets.HBox([count, coupling, spin, remove])
        remove.on_click(lambda _: setattr(self.rows, "children", tuple(child for child in self.rows.children if child is not row)))
        self.rows.children = (*self.rows.children, row)

    def set(self, radical):
        self.g.value = radical.g
        self.rows.children = ()
        for group in radical.groups:
            self.append(group)

    def value(self):
        return Radical(self.g.value, [Group(*(control.value for control in row.children[:3])) for row in self.rows.children])


class Workbench:
    def __init__(self):
        self.config = Config()
        self.shape = None
        self.source = "simulation"
        self.metadata = {}
        self.result = simulate(self.config)
        self.cursor_index = 0
        self.controls = {}
        self.shape_controls = {}
        self.editors = [RadicalEditor(self.config.radical1), RadicalEditor(self.config.radical2)]
        radicals = widgets.Tab(children=[editor.widget for editor in self.editors])
        for index in range(2):
            radicals.set_title(index, f"Radical {index + 1}")
        config_fields = [self.control("identical", "Identical radicals", True),
                         self.control("second_order", "Second-order splitting", True),
                         self.control("centre_resonance", "Centre on resonance", True),
                         self.control("centre", "Centre / G", 3483.0),
                         self.control("sweep", "Sweep / G", 180.0),
                         self.control("frequency", "Frequency / GHz", 9.777),
                         self.control("points", "Points", 1024),
                         self.control("pair", "Pair", "T", ["S", "T", "F"]),
                         self.control("observed_radical", "Observed radical", 1, [1, 2]),
                         self.control("linear", "Linear coefficient", 0.0),
                         self.control("use_exponent", "Use Q exponent", False),
                         self.control("exponent", "Q exponent", 0.5),
                         self.control("net_effect", "ST-minus net effect", True),
                         self.control("normalize", "Normalize components", True)]
        mixing = [self.control("rpm_weight", "ST-zero RPM", 1.0),
                  self.control("stm_weight", "ST-minus RPM", 0.0),
                  self.control("tm_weight", "Triplet mechanism", 0.0),
                  self.control("tm_ratio", "TM polarization ratio", -50.0)]
        self.mode = widgets.Dropdown(options=["Sticks", "Relaxation", "Slow exchange"], description="Calculation")
        shaping = [self.mode]
        for name, label in [("t1", "T1 / us"), ("t2", "T2 / us"), ("microwave", "Microwave / rad MHz"),
                            ("delay", "Delay / us"), ("interval", "Interval / us"), ("samples", "Samples"),
                            ("points", "Lineshape points"), ("exchange_time", "Exchange time / us"),
                            ("polarization", "Exchange polarization")]:
            value = getattr(Shape(), name)
            shaping.append(self.control(name, label, 5.5 if value is None else value, target=self.shape_controls))
        sections = widgets.Accordion(children=[radicals, widgets.VBox(config_fields), widgets.VBox(mixing), widgets.VBox(shaping)])
        sections.selected_index = 0
        for index, title in enumerate(["Radicals", "Field and physics", "Mechanism mixture", "Lineshape"]):
            sections.set_title(index, title)
        self.run_button = widgets.Button(description="Calculate", icon="play", button_style="primary")
        self.run_button.on_click(lambda _: asyncio.get_running_loop().create_task(self.calculate()))
        reset = widgets.Button(icon="undo", tooltip="Restore defaults", layout=widgets.Layout(width="36px"))
        reset.on_click(self.reset)
        self.status = widgets.HTML()
        self.filename = widgets.Text(value="spectrum.json", description="JSON file")
        save_button = widgets.Button(icon="save", tooltip="Save JSON", layout=widgets.Layout(width="36px"))
        save_button.on_click(self.save)
        self.upload = widgets.FileUpload(accept=".json,.bin,.dat", multiple=False, description="Open / Import")
        self.upload.observe(self.open, names="value")
        self.view = widgets.Dropdown(options=["mixed", "rpm", "stm", "tm", "equilibrium"], value="mixed", description="Spectrum")
        self.view.observe(lambda _: self.draw(), names="value")
        self.readout = widgets.HTML()
        with plt.ioff():
            self.figure, self.axes = plt.subplots(figsize=(6, 4.2), layout="constrained")
        if isinstance(self.figure.canvas, widgets.Widget):
            self.figure.canvas.layout.width = "100%"
            self.figure.canvas.layout.min_width = "0px"
            self.figure.canvas.header_visible = False
            self.figure.canvas.toolbar_position = "bottom"
        self.figure.canvas.mpl_connect("button_press_event", self.click)
        self.figure.canvas.mpl_connect("key_press_event", self.key)
        self.plot_output = widgets.Output(layout=widgets.Layout(width="100%", min_width="0px", overflow="auto"))
        with self.plot_output:
            display(self.figure.canvas if isinstance(self.figure.canvas, widgets.Widget) else self.figure)
        self.widget = widgets.VBox([
            widgets.HTML("<h2 style='font-family:Georgia;color:#12695d;margin:8px 0'>CIDEP</h2>"),
            widgets.HBox([self.run_button, reset, self.status], layout=widgets.Layout(flex_flow="row wrap")),
            widgets.HBox([widgets.VBox([sections], layout=widgets.Layout(width="350px", max_width="100%")),
                          widgets.VBox([self.view, self.plot_output, self.readout], layout=widgets.Layout(flex="1 1 350px", min_width="0px", max_width="100%"))],
                         layout=widgets.Layout(flex_flow="row wrap", width="100%")),
            widgets.HBox([self.upload, self.filename, save_button], layout=widgets.Layout(flex_flow="row wrap"))])
        self.draw()

    def control(self, name, label, value, options=None, target=None):
        target = self.controls if target is None else target
        arguments = dict(value=value, description=label, style={"description_width": "165px"},
                         layout=widgets.Layout(width="335px", max_width="100%"))
        if options is not None:
            control = widgets.Dropdown(options=options, **arguments)
        elif isinstance(value, bool):
            control = widgets.Checkbox(**arguments)
        elif isinstance(value, int):
            control = widgets.IntText(**arguments)
        else:
            control = widgets.FloatText(**arguments)
        target[name] = control
        return control

    def inputs(self):
        values = {key: control.value for key, control in self.controls.items()}
        if not values.pop("use_exponent"):
            values["exponent"] = None
        config = Config(radical1=self.editors[0].value(), radical2=self.editors[1].value(), **values)
        config.validate()
        shape = None
        if self.mode.value != "Sticks":
            values = {key: control.value for key, control in self.shape_controls.items()}
            if self.mode.value != "Slow exchange":
                values["exchange_time"] = None
            shape = Shape(**values)
            shape.validate()
        return config, shape

    async def calculate(self):
        self.run_button.disabled = True
        self.status.value = "Calculating..."
        try:
            config, shape = self.inputs()
            result = await asyncio.to_thread(simulate, config)
            if shape:
                result = await asyncio.to_thread(broaden, result, shape)
            self.config, self.shape, self.result, self.source = config, shape, result, "simulation"
            self.metadata = {}
            self.refresh_views()
            self.status.value = f"{config.points} points"
        except Exception as error:
            self.status.value = html.escape(str(error))
        finally:
            self.run_button.disabled = False

    def refresh_views(self):
        self.view.options = [key for key in self.result if key in ARRAYS and key != "field"]
        self.view.value = "shaped" if "shaped" in self.result else "mixed"
        self.draw()

    def draw(self):
        if not hasattr(self, "axes") or self.view.value not in self.result:
            return
        self.axes.clear()
        fields, values = self.result["field"], self.result[self.view.value]
        self.axes.axhline(0, color="#c9cecb", linewidth=0.7)
        self.axes.plot(fields, values, color="#087f72", linewidth=1)
        self.axes.set(xlabel="Magnetic field / Gauss", ylabel="Signal", title=self.view.value.replace("_", " ").title())
        self.axes.grid(alpha=0.15)
        self.cursor_index = min(self.cursor_index, len(fields) - 1)
        self.cursor = self.axes.axvline(fields[self.cursor_index], color="#c43f52", linewidth=0.8)
        self.update_cursor()

    def update_cursor(self):
        fields, values = self.result["field"], self.result[self.view.value]
        self.cursor.set_xdata([fields[self.cursor_index]] * 2)
        self.readout.value = f"<code>{self.cursor_index + 1} / {len(fields)} &nbsp; {fields[self.cursor_index]:.6f} G &nbsp; {values[self.cursor_index]:.7g}</code>"
        self.figure.canvas.draw_idle()

    def click(self, event):
        if event.inaxes is self.axes and event.xdata is not None:
            self.cursor_index = int(np.argmin(np.abs(self.result["field"] - event.xdata)))
            self.update_cursor()

    def key(self, event):
        if event.key in ("left", "right"):
            self.cursor_index = max(0, min(len(self.result["field"]) - 1, self.cursor_index + (1 if event.key == "right" else -1)))
            self.update_cursor()

    def set_config(self, config, shape=None):
        for editor, radical in zip(self.editors, [config.radical1, config.radical2]):
            editor.set(radical)
        for key, value in asdict(config).items():
            if key in self.controls and value is not None:
                self.controls[key].value = value
        self.controls["use_exponent"].value = config.exponent is not None
        self.mode.value = "Sticks" if shape is None else "Slow exchange" if shape.exchange_time is not None else "Relaxation"
        for key, value in asdict(shape or Shape()).items():
            if value is not None:
                self.shape_controls[key].value = value

    def reset(self, _):
        self.set_config(Config())
        self.status.value = "Defaults restored"

    def save(self, _):
        try:
            save(self.filename.value, self.result, self.config, self.shape, self.source, self.metadata)
            self.status.value = "JSON saved"
        except Exception as error:
            self.status.value = html.escape(str(error))

    def open(self, change):
        if not change["new"]:
            return
        try:
            uploaded = change["new"][0]
            raw = bytes(uploaded["content"])
            payload = json.loads(raw.decode("utf-8")) if uploaded["name"].lower().endswith(".json") else import_legacy(raw)
            validate(payload)
            self.result = {key: np.asarray(value, dtype=float) if key in ARRAYS else value for key, value in payload["data"].items()}
            self.config = Config.from_dict(payload["config"]) if payload.get("config") else None
            self.shape = Shape(**payload["shape"]) if payload.get("shape") else None
            self.source = payload.get("source", "import")
            self.metadata = payload.get("metadata", {})
            if self.config:
                self.set_config(self.config, self.shape)
            self.refresh_views()
            self.status.value = "Spectrum loaded"
        except Exception as error:
            self.status.value = html.escape(str(error))


def launch():
    workbench = Workbench()
    display(workbench.widget)
    return workbench