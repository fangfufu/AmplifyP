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
        
        # Wait for semantics
        page.wait_for_selector("flt-semantics-host", state="attached", timeout=30000)
        placeholder = page.locator("flt-semantics-placeholder").first
        placeholder.dispatch_event("click")
        
        # Wait for inputs
        page.wait_for_selector('input:not([type="file"])', state="attached", timeout=30000)
        
        # V1
        inputs = page.locator('input:not([type="file"])')
        inputs.nth(0).click(force=True)
        inputs.nth(0).press_sequentially("V1")
        page.keyboard.press("Tab")
        inputs.nth(1).press_sequentially("ATGCATGCATGCATGC")
        page.keyboard.press("Enter")
        time.sleep(1)
        
        # V2
        inputs = page.locator('input:not([type="file"])')
        inputs.nth(2).click(force=True)
        inputs.nth(2).press_sequentially("V2")
        page.keyboard.press("Tab")
        inputs.nth(3).press_sequentially("GCATGCATGCATGCAT")
        page.keyboard.press("Enter")
        time.sleep(1)
        
        # I1
        inputs = page.locator('input:not([type="file"])')
        inputs.nth(4).click(force=True)
        inputs.nth(4).press_sequentially("I1")
        page.keyboard.press("Tab")
        inputs.nth(5).press_sequentially("XYZXYZXYZXYZ")
        page.keyboard.press("Enter")
        time.sleep(2)
        
        # Click last input
        page.locator('input:not([type="file"])').last.click(force=True)
        time.sleep(1)
        
        # Print initial input count
        initial_count = page.locator('input:not([type="file"])').count()
        print("Initial count before manual add click:", initial_count)
        
        # Find the Add button using role search inside the merged row element
        row = page.locator("[aria-label*='Add Primer Below']").first
        add_btn = row.locator("[role='button']").first
        
        print("Add button visible:", add_btn.is_visible())
        print("Add button enabled:", add_btn.is_enabled())
        print("Add button bounding box:", add_btn.bounding_box())
        
        # Click it
        add_btn.click(force=True)
        time.sleep(2)
        
        # Print final count
        final_count = page.locator('input:not([type="file"])').count()
        print("Final count after manual add click:", final_count)
        
        browser.close()
    server.shutdown()
    thread.join()

if __name__ == "__main__":
    main()
