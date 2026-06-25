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
        
        page.goto(f"http://localhost:{SERVER_PORT}/?enable-semantics=true")
        
        # Wait for page title
        page.wait_for_function("document.title === 'AmplifyP'", timeout=60000)
        
        # Wait 20 seconds for Python worker to completely initialize and UI to render
        time.sleep(20)
        
        # Check if placeholder is there and click it
        placeholder = page.locator("flt-semantics-placeholder").first
        if placeholder.count() > 0:
            placeholder.dispatch_event("click")
            
        time.sleep(5)
        
        # Print all textareas
        textareas = page.locator('textarea')
        print(f"Total textarea elements found: {textareas.count()}")
        for i in range(textareas.count()):
            label = textareas.nth(i).evaluate('el => el.getAttribute("aria-label")')
            outer = textareas.nth(i).evaluate('el => el.outerHTML')
            print(f"Textarea {i}: label = '{label}', outerHTML = {outer}")
            
        browser.close()
    server.shutdown()
    thread.join()

if __name__ == "__main__":
    main()
