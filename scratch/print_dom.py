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
        
        # Wait a bit for loading
        time.sleep(5)
        print("Page Title:", page.title())
        print("HTML content of body:")
        body_content = page.eval_on_selector("body", "el => el.innerHTML")
        print(body_content)
        
        browser.close()
    server.shutdown()
    thread.join()

if __name__ == "__main__":
    main()
