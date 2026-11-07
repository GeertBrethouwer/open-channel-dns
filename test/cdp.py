"""Evaluate JavaScript in the Android app's WebView over the Chrome DevTools Protocol.

Debug builds expose a DevTools socket. Forward it first:
    adb forward tcp:9333 localabstract:webview_devtools_remote_$(adb shell pidof io.github.geertbrethouwer.openchanneldns)
Then:
    python3 test/cdp.py 9333 'document.title'
No third-party packages needed.
"""
import base64, json, os, socket, struct, sys, urllib.request

port = int(sys.argv[1])
expr = sys.argv[2]
pages = json.load(urllib.request.urlopen(f"http://127.0.0.1:{port}/json"))
ws_url = next(p["webSocketDebuggerUrl"] for p in pages if p.get("type") == "page")
path = "/" + ws_url.split("/", 3)[3]

s = socket.create_connection(("127.0.0.1", port), timeout=120)
key = base64.b64encode(os.urandom(16)).decode()
s.sendall((f"GET {path} HTTP/1.1\r\nHost: 127.0.0.1:{port}\r\nUpgrade: websocket\r\nConnection: Upgrade\r\n"
           f"Sec-WebSocket-Key: {key}\r\nSec-WebSocket-Version: 13\r\n\r\n").encode())
resp = b""
while b"\r\n\r\n" not in resp:
    resp += s.recv(4096)


def send(obj):
    data = json.dumps(obj).encode()
    mask = os.urandom(4)
    hdr = bytearray([0x81])
    n = len(data)
    if n < 126:
        hdr.append(0x80 | n)
    elif n < 65536:
        hdr.append(0x80 | 126); hdr += struct.pack(">H", n)
    else:
        hdr.append(0x80 | 127); hdr += struct.pack(">Q", n)
    s.sendall(bytes(hdr) + mask + bytes(b ^ mask[i % 4] for i, b in enumerate(data)))


def recv_exact(n):
    buf = b""
    while len(buf) < n:
        chunk = s.recv(n - len(buf))
        if not chunk:
            raise EOFError
        buf += chunk
    return buf


def recv():
    msg = b""
    while True:
        b1, b2 = recv_exact(2)
        n = b2 & 0x7F
        if n == 126: n = struct.unpack(">H", recv_exact(2))[0]
        elif n == 127: n = struct.unpack(">Q", recv_exact(8))[0]
        msg += recv_exact(n)
        if b1 & 0x80:
            return json.loads(msg)


send({"id": 1, "method": "Runtime.evaluate",
      "params": {"expression": expr, "awaitPromise": True, "returnByValue": True}})
while True:
    m = recv()
    if m.get("id") == 1:
        r = m.get("result", {})
        if "exceptionDetails" in r:
            print("EXCEPTION:", json.dumps(r["exceptionDetails"])[:1500])
        else:
            v = r.get("result", {}).get("value")
            print(v if isinstance(v, str) else json.dumps(v))
        break
