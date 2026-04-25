package org.cenicana.bio.web;

import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;
import org.cenicana.bio.cli.VcfToolkit;

import java.io.*;
import java.net.InetSocketAddress;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Map;
import java.util.Scanner;

public class BioCenicanaWebServer {

    public static void start(int port) throws IOException {
        HttpServer server = HttpServer.create(new InetSocketAddress(port), 0);
        
        // Static files handler (index.html)
        server.createContext("/", new HttpHandler() {
            @Override
            public void handle(HttpExchange exchange) throws IOException {
                String path = exchange.getRequestURI().getPath();
                if (path.equals("/")) path = "/index.html";
                
                // 1. Try static resources first
                String resourcePath = "src/main/resources/static" + path;
                File file = new File(resourcePath);
                
                // 2. If not found, try project root (for results)
                if (!file.exists()) {
                    file = new File("." + path);
                }
                
                if (file.exists() && !file.isDirectory()) {
                    byte[] response = Files.readAllBytes(file.toPath());
                    String contentType = "text/plain";
                    if (path.endsWith(".html")) contentType = "text/html";
                    if (path.endsWith(".css")) contentType = "text/css";
                    if (path.endsWith(".js")) contentType = "application/javascript";
                    if (path.endsWith(".png")) contentType = "image/png";
                    
                    exchange.getResponseHeaders().set("Content-Type", contentType);
                    exchange.sendResponseHeaders(200, response.length);
                    OutputStream os = exchange.getResponseBody();
                    os.write(response);
                    os.close();
                } else {
                    String error = "404 Not Found";
                    exchange.sendResponseHeaders(404, error.length());
                    exchange.getResponseBody().write(error.getBytes());
                    exchange.getResponseBody().close();
                }
            }
        });

        // API handler
        server.createContext("/api/analysis/run", new HttpHandler() {
            @Override
            public void handle(HttpExchange exchange) throws IOException {
                if (!exchange.getRequestMethod().equalsIgnoreCase("POST")) {
                    exchange.sendResponseHeaders(405, -1);
                    return;
                }

                InputStream is = exchange.getRequestBody();
                String body = new Scanner(is, StandardCharsets.UTF_8).useDelimiter("\\A").next();
                
                // Simple JSON parsing (since we have no Jackson/Gson)
                String command = extractValue(body, "command");
                String argsStr = extractValue(body, "args");
                
                String[] args = argsStr.split("\\s+");
                String[] fullArgs = new String[args.length + 1];
                fullArgs[0] = command;
                System.arraycopy(args, 0, fullArgs, 1, args.length);

                ByteArrayOutputStream baos = new ByteArrayOutputStream();
                PrintStream ps = new PrintStream(baos);
                PrintStream oldOut = System.out;
                
                String responseJson;
                try {
                    System.setOut(ps);
                    VcfToolkit.main(fullArgs);
                    System.out.flush();
                    responseJson = "{\"success\": true, \"output\": \"" + 
                        baos.toString().replace("\"", "\\\"").replace("\n", "\\n").replace("\r", "") + "\"}";
                } catch (Exception e) {
                    responseJson = "{\"success\": false, \"error\": \"" + e.getMessage() + "\"}";
                } finally {
                    System.setOut(oldOut);
                }

                exchange.getResponseHeaders().set("Content-Type", "application/json");
                exchange.sendResponseHeaders(200, responseJson.getBytes().length);
                OutputStream os = exchange.getResponseBody();
                os.write(responseJson.getBytes());
                os.close();
            }
        });

        System.out.println("[WEB] BioCenicana Web Server started at http://localhost:" + port);
        server.setExecutor(null); 
        server.start();
    }

    private static String extractValue(String json, String key) {
        int start = json.indexOf("\"" + key + "\":") + key.length() + 3;
        int end = json.indexOf("\"", start + 1);
        if (json.charAt(start) == '\"') {
            return json.substring(start + 1, end);
        } else {
            // Numeric or unquoted value
            int endAlt = json.indexOf(",", start);
            if (endAlt == -1) endAlt = json.indexOf("}", start);
            return json.substring(start, endAlt).trim();
        }
    }
}
