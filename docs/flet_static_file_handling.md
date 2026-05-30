# Static Flet Deployment File Handling

This document explains the mechanism used to enable file uploads and downloads
when deploying a Flet application as a static website (Pyodide).

## The Modern Solution: Flet Native FilePicker

Previously, Flet required custom hot-patching workarounds to bridge files
between the Pyodide Web Worker environment and the browser's main thread. With
modern Flet updates, this bridging is handled completely **natively** under the
hood by Flet's standard `FilePicker` API.

### 1. File Upload (Opening / Loading Files)

To read a file in a static or web deployment without a traditional server
backend, we utilize `pick_files()` with `with_data=True`. This loads the
selected file's content directly into memory in the browser.

#### implementation in `app.py`

```python
files = await file_picker.pick_files(
    dialog_title="Load State",
    allowed_extensions=["yaml", "yml"],
    file_type=ft.FilePickerFileType.CUSTOM,
    with_data=True,
)
if not files:
    return

file = files[0]
if file.bytes:
    # Access file bytes directly from memory (supported on web/Pyodide)
    content = file.bytes.decode("utf-8")
    parsed_state = yaml.safe_load(content)
else:
    # Fallback to reading file path (supported on desktop)
    with open(file.path) as f:
        parsed_state = yaml.safe_load(f)
```

Adding `with_data=True` triggers Flet to read the file's data into a `bytes`
attribute on the returned `FilePickerFile` object, bypassing the lack of local
file system path access in the browser sandbox.

### 2. File Download (Saving Files)

To download generated state/data as a file on web/static platforms, we provide
the file's content directly to the `save_file()` dialog using the `src_bytes`
parameter.

#### Implementation in `app.py`

```python
file_path = await file_picker.save_file(
    dialog_title="Save State",
    file_name=STATE_FILE,
    allowed_extensions=["yaml", "yml"],
    file_type=ft.FilePickerFileType.CUSTOM,
    src_bytes=yaml_str.encode("utf-8"),
)
if page.web:
    # Flet Web handles generating a Blob and downloading it natively
    show_snackbar("State ready for download!")
else:
    # On desktop, save_file only returns the path, and we write to it
    if file_path is None:
        return
    with open(file_path, "w") as f:
        f.write(yaml_str)
    show_snackbar("State saved successfully!")
```

## Deployment Workflow

Because file picking and saving are now native to the Flet framework, our
deployment workflow is incredibly streamlined:

1. **Build:** `flet publish src/main.py` compiles the application into static
   HTML/JS/WASM assets in `src/dist`.
1. **Deploy:** The assets are deployed directly to a static server (e.g., GitHub
   Pages) without any hot-patching scripts.
