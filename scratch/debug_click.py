import os
import time
from http.server import HTTPServer, SimpleHTTPRequestHandler
import threading
from playwright.sync_api import sync_playwright

SERVER_PORT = 34521
DIST_DIR = os.path.join(os.getcwd(), "dist")

def serve_app():
    HTTPServer.allow_reuse_address = True
    class CustomHTTPRequestHandler(SimpleHTTPRequestHandler):
        def translate_path(self, path: str) -> str:
            if path.startswith("/AmplifyP/"):
                path = path[len("/AmplifyP") :]
            elif path == "/AmplifyP":
                path = "/"
            return super().translate_path(path)

    server = HTTPServer(
        ("localhost", SERVER_PORT),
        lambda *args: CustomHTTPRequestHandler(*args, directory=DIST_DIR),
    )
    thread = threading.Thread(target=server.serve_forever)
    thread.daemon = True
    thread.start()
    return server, thread

def main():
    server, thread = serve_app()
    time.sleep(1)
    
    with sync_playwright() as p:
        browser = p.chromium.launch(
            headless=True,
            args=[
                "--disable-gpu",
                "--no-sandbox",
                "--disable-dev-shm-usage",
                "--ignore-gpu-blocklist",
            ]
        )
        page = browser.new_page()
        page.set_default_timeout(30000)
        page.on("console", lambda msg: print(f"Browser console: {msg.text}"))
        
        print("Navigating to page...")
        page.goto(f"http://localhost:{SERVER_PORT}/?enable-semantics=true")
        
        # Wait for page title
        print("Waiting for title 'AmplifyP'...")
        page.wait_for_function("document.title === 'AmplifyP'", timeout=60000)
        print("Title acquired:", page.title())
        
        # Let it settle for a couple seconds
        time.sleep(3)
        
        # Print semantics host state
        print("flt-semantics-host content before click:")
        host_html = page.eval_on_selector("flt-semantics-host", "el => el.innerHTML")
        print(host_html)
        
        # Click the placeholder
        placeholder = page.locator("flt-semantics-placeholder").first
        print("Clicking placeholder using dispatch_event...")
        placeholder.dispatch_event("click")
        time.sleep(1)
        print("flt-semantics-host content after dispatch_event click:")
        host_html = page.eval_on_selector("flt-semantics-host", "el => el.innerHTML")
        print(host_html)
        
        print("Clicking placeholder using click(force=True)...")
        placeholder.click(force=True, timeout=5000)
        time.sleep(2)
        print("flt-semantics-host content after click(force=True):")
        host_html = page.eval_on_selector("flt-semantics-host", "el => el.innerHTML")
        print(host_html)
        
        browser.close()
    server.shutdown()
    thread.join()

if __name__ == "__main__":
    main()
