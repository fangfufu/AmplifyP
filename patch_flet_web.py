"""Patches Flet static web files to support FilePicker in Pyodide/Web."""

import os
import sys


def patch_flet_web(dist_dir: str) -> None:
    """Patch python.js and python-worker.js in the distribution directory."""
    python_js_path = os.path.join(dist_dir, "python.js")
    python_worker_path = os.path.join(dist_dir, "python-worker.js")

    if not os.path.exists(python_js_path) or not os.path.exists(
        python_worker_path
    ):
        print("Could not find Flet static web files to patch.")
        sys.exit(1)

    # Patch python.js
    with open(python_js_path, encoding="utf-8") as f:
        js_code = f.read()

    # JS code is indented to match context, but Python string literal doesn't
    # enforce indentation. We can keep it simple.
    # Note: Long lines in JS string are unavoidable unless we break them up,
    # but ruff checks Python line length.
    interceptor = (
        "\n"
        "        if (typeof event.data === 'string') {\n"
        "            if (event.data === 'open_file_picker') {\n"
        "                console.log('[PATCH] open_file_picker triggered');\n"
        "                let input = document.createElement('input');\n"
        "                input.type = 'file';\n"
        "                input.accept = '.yaml,.yml';\n"
        "                input.style.display = 'none';\n"
        "                input.onchange = (e) => {\n"
        "                    let file = e.target.files[0];\n"
        "                    if (!file) return;\n"
        "                    let reader = new FileReader();\n"
        "                    reader.onload = (e2) => {\n"
        "                        app.worker.postMessage("
        "{type: 'file_content', content: e2.target.result});\n"
        "                    };\n"
        "                    reader.readAsText(file);\n"
        "                };\n"
        "                input.click();\n"
        "                return;\n"
        "            }\n"
        "            try {\n"
        "                let msg = JSON.parse(event.data);\n"
        "                if (msg && msg.type === 'save_file') {\n"
        "                    console.log('[PATCH] save_file triggered', "
        "msg.filename);\n"
        "                    let blob = new Blob([msg.content], "
        "{ type: 'text/yaml' });\n"
        "                    let url = URL.createObjectURL(blob);\n"
        "                    let a = document.createElement('a');\n"
        "                    a.href = url;\n"
        "                    a.download = msg.filename;\n"
        "                    a.style.display = 'none';\n"
        "                    document.body.appendChild(a);\n"
        "                    a.click();\n"
        "                    document.body.removeChild(a);\n"
        "                    setTimeout(() => "
        "URL.revokeObjectURL(url), 1000);\n"
        "                    return;\n"
        "                }\n"
        "            } catch (e) {}\n"
        "        }\n"
    )

    if 'event.data === "open_file_picker"' not in js_code:
        # Insert BEFORE the string check
        js_code = js_code.replace(
            'if (typeof event.data === "string") {',
            interceptor + 'if (typeof event.data === "string") {',
        )
        with open(python_js_path, "w", encoding="utf-8") as f:
            f.write(js_code)
        print("Patched python.js for FilePicker support.")
    else:
        print("python.js already patched.")

    # Patch python-worker.js
    with open(python_worker_path, encoding="utf-8") as f:
        worker_code = f.read()

    worker_patch = (
        "\n"
        "        if (event.data && event.data.type === 'file_content') {\n"
        "            if (self.custom_file_callback) "
        "self.custom_file_callback(event.data.content);\n"
        "        } else {\n"
        "            flet_js.send(event.data);\n"
        "        }\n"
    )

    if "event.data.type === 'file_content'" not in worker_code:
        worker_code = worker_code.replace(
            "flet_js.send(event.data);", worker_patch
        )
        with open(python_worker_path, "w", encoding="utf-8") as f:
            f.write(worker_code)
        print("Patched python-worker.js for FilePicker support.")
    else:
        print("python-worker.js already patched.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python patch_flet_web.py <dist_dir>")
        sys.exit(1)
    patch_flet_web(sys.argv[1])
