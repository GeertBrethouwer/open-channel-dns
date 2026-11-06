package io.github.geertbrethouwer.openchanneldns;

import android.annotation.SuppressLint;
import android.app.Activity;
import android.os.Bundle;
import android.view.WindowManager;
import android.webkit.WebSettings;
import android.webkit.WebView;

/** Shows the same page as the website in a WebView, with the C++ solver behind it. */
public class MainActivity extends Activity {
    private WebView web;
    private final SolverBridge bridge = new SolverBridge();

    @SuppressLint({"SetJavaScriptEnabled", "AddJavascriptInterface"})
    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        getWindow().addFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON);
        web = new WebView(this);
        WebSettings s = web.getSettings();
        s.setJavaScriptEnabled(true);
        s.setDomStorageEnabled(true);
        s.setAllowFileAccess(false);
        s.setAllowContentAccess(false);
        // The page detects window.OCDNSNative and switches from the JavaScript solver to C++.
        web.addJavascriptInterface(bridge, "OCDNSNative");
        setContentView(web);
        if (savedInstanceState != null) web.restoreState(savedInstanceState);
        else web.loadUrl("file:///android_asset/index.html");
    }

    @Override
    protected void onSaveInstanceState(Bundle out) {
        super.onSaveInstanceState(out);
        web.saveState(out);
    }

    @Override
    protected void onPause() {
        super.onPause();
        bridge.setSuspended(true);
        web.onPause();
    }

    @Override
    protected void onResume() {
        super.onResume();
        web.onResume();
        bridge.setSuspended(false);
    }

    @Override
    protected void onDestroy() {
        web.destroy();
        bridge.destroy();
        super.onDestroy();
    }
}
