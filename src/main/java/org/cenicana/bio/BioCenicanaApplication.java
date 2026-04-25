package org.cenicana.bio;

import org.cenicana.bio.cli.VcfToolkit;
import org.cenicana.bio.web.BioCenicanaWebServer;
import java.io.IOException;

public class BioCenicanaApplication {

    public static void main(String[] args) {
        if (args.length > 0 && !args[0].equals("web")) {
            // Standard CLI Mode
            VcfToolkit.main(args);
        } else {
            // Web Server Mode
            try {
                BioCenicanaWebServer.start(9090);
            } catch (IOException e) {
                System.err.println("[ERROR] Could not start web server: " + e.getMessage());
            }
        }
    }
}
