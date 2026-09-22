"use strict";

// get all of the releases from versions.json, and use these to populate the
// dropdown menu of different releases
$(document).ready(function () {
    // Define versions_json_url in versions.html through a template variable
    $.getJSON(versions_json_url)
        .done(function (data) {
            
            // Try to determine the path of the current page relative to the doc root
            var page_path = "";
            var current_url = window.location.href;
            
            // Sort data by URL length descending to match the most specific first
            var sorted_by_length = data.slice().sort(function (a, b) {
                return b.url.length - a.url.length;
            });
            
            // Check if the current URL starts with one of the known version URLs
            $.each(sorted_by_length, function (i, item) {
                if (current_url.startsWith(item.url)) {
                    page_path = current_url.substring(item.url.length);
                    return false; // break
                }
            });
            
            // Fallback to DOCUMENTATION_OPTIONS if available (e.g. for local browsing)
            if (!page_path && typeof DOCUMENTATION_OPTIONS !== 'undefined' && DOCUMENTATION_OPTIONS.pagename) {
                var pagename = DOCUMENTATION_OPTIONS.pagename;
                page_path = "/" + (pagename === "index" ? "" : pagename + ".html") + window.location.search + window.location.hash;
            }

            // Ensure page_path starts with a slash
            if (page_path && !page_path.startsWith("/")) {
                page_path = "/" + page_path;
            }

            // Populate the dropdown
            $.each(data.sort(function (a, b) {
                return a.version > b.version;
            }), function (i, item) {
                var target_url = item.url;
                if (target_url.endsWith("/")) {
                    target_url = target_url.slice(0, -1);
                }
                target_url += page_path;

                $("<dd>").append(
                    $("<a>").text(item.display).attr('href', target_url)
                ).appendTo("#versionselector");
            });
        })
        .fail(function (d, textStatus, error) {
            console.error("getJSON failed, status: " + textStatus + ", error: " + error);
        });
});

console.log("Loading versions from " + versions_json_url);
