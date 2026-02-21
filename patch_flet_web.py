import os
import sys

def patch_flet_web(dist_dir):
    python_js_path = os.path.join(dist_dir, "python.js")
    python_worker_path = os.path.join(dist_dir, "python-worker.js")

    if not os.path.exists(python_js_path) or not os.path.exists(python_worker_path):
        print("Could not find Flet static web files to patch.")
        sys.exit(1)

    # Patch python.js
    with open(python_js_path, "r", encoding="utf-8") as f:
        js_code = f.read()

    interceptor = """
        if (typeof event.data === 'string' && event.data === 'open_file_picker') {
            let input = document.createElement('input');
            input.type = 'file';
            input.accept = '.yaml,.yml';
            input.onchange = (e) => {
                let file = e.target.files[0];
                if (!file) return;
                let reader = new FileReader();
                reader.onload = (e2) => {
                    app.worker.postMessage({type: 'file_content', content: e2.target.result});
                };
                reader.readAsText(file);
            };
            input.click();
            return;
        }
        if (typeof event.data === 'object' && event.data.type === 'save_file') {
            let blob = new Blob([event.data.content], { type: 'text/yaml' });
            let url = URL.createObjectURL(blob);
            let a = document.createElement('a');
            a.href = url;
            a.download = event.data.filename;
            a.click();
            URL.revokeObjectURL(url);
            return;
        }
"""
    if 'event.data === "open_file_picker"' not in js_code:
        js_code = js_code.replace('if (typeof event.data === "string") {', 'if (typeof event.data === "string") {' + interceptor)
        with open(python_js_path, "w", encoding="utf-8") as f:
            f.write(js_code)
        print("Patched python.js for FilePicker support.")
    else:
        print("python.js already patched.")

    # Patch python-worker.js
    with open(python_worker_path, "r", encoding="utf-8") as f:
        worker_code = f.read()

    worker_patch = """
        if (event.data && event.data.type === 'file_content') {
            if (self.custom_file_callback) self.custom_file_callback(event.data.content);
        } else {
            flet_js.send(event.data);
        }
"""
    if 'event.data.type === \'file_content\'' not in worker_code:
        worker_code = worker_code.replace("flet_js.send(event.data);", worker_patch)
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
