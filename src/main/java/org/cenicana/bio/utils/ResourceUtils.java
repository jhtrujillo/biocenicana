package org.cenicana.bio.utils;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;

/**
 * Utility to load resources from the classpath.
 */
public class ResourceUtils {

    /**
     * Loads a resource as a String.
     * @param resourceName The name of the resource (e.g., "plotly.min.js")
     * @return The resource content as a UTF-8 string, or an empty string if not found.
     */
    public static String loadResource(String resourceName) {
        try (InputStream is = ResourceUtils.class.getClassLoader().getResourceAsStream(resourceName)) {
            if (is == null) {
                System.err.println("  [Warning] Resource not found: " + resourceName);
                return "";
            }
            return new String(is.readAllBytes(), StandardCharsets.UTF_8);
        } catch (IOException e) {
            System.err.println("  [Error] Failed to load resource " + resourceName + ": " + e.getMessage());
            return "";
        }
    }
}
