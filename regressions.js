'use strict';

$(document).ready(function() {
    /* Cached contents of downloaded regressions.json */
    var regression_data = null;
    /* Current page title */
    var current_title = "All regressions";
    /* Whether HTML5 local storage is available */
    var local_storage_available = false;
    /* Key prefix for ignored regressions. For each ignored regression,
       a key "ignore_key_prefix + md5(benchmark_name + date_a + date_b)"
       is added to HTML5 local storage.
     */
    var ignore_key_prefix = null;
    /* Set of ignored regressions, same information as in HTML5 local storage.
       Useful if local storage runs out of space. */
    var ignored_regressions = {};

    function load_data(params) {
        $("#title").text(current_title);

        if (typeof(Storage) !== "undefined") {
            /* html5 local storage available */
            local_storage_available = true;
        }

        if (regression_data !== null) {
            // already displayed
        }
        else {
            var message = $('<div>Loading...</div>');
            $('#regressions-body').append(message);
            $.ajax({
                url: 'regressions.json',
                cache: false
            }).done(function (data) {
                regression_data = data;
                var main_div = display_data(data, params);
                $('#regressions-body').empty();
                $('#regressions-body').append(main_div);
            });
        }
    }

    function display_data(data, params) {
        var main_div = $('<div/>');
        var branches = $.asv.master_json.params['branch'];
        var all_ignored_keys = {};

        ignore_key_prefix = 'asv-r-' + $.asv.master_json.project;

        if (branches && branches.length > 1) {
            /* Add a branch selector */
            var dropdown_menu = $('<ul class="dropdown-menu" role="menu"/>');
            var dropdown_div = $('<div class="dropdown">');

            dropdown_div.append($('<button class="btn btn-default dropdown-toggle" data-toggle="dropdown">Branches ' +
                                  '<span class="caret"/></button>'));
            dropdown_div.append(dropdown_menu);

            main_div.append(dropdown_div);
        }

        var feed_div = $('<div class="feed-div"><a class="btn" href="regressions.xml">Feed (Atom)</a></div>');
        main_div.append(feed_div);

        $.each(branches, function(i, branch) {
            var branch_div = $('<div class="regression-div"/>')

            var display_table = $('<table class="table table-hover"/>');
            var ignored_table = $('<table class="table table-hover ignored"/>');
            var ignored_button = $('<button class="btn btn-default">Show ignored regressions...</button>');
            var ignored_conf_sample_div = $('<div class="ignored"/>');

            if (branches && branches.length > 1) {
                var branch_link = $('<a/>')
                branch_link.text(branch);

                dropdown_menu.append($('<li role="presentation"/>').append(branch_link));
                branch_link.on('click', function(evt) {
                    current_title = "Regressions in " + branch + " branch";
                    $("#title").text(current_title);
                    $(".regression-div").hide();
                    $(".ignored").hide();
                    ignored_button.show();
                    $("#regression-div-" + i).show();
                    $("#regression-div-" + i + '-ignored').show();
                });
            }
            else {
                branch = null;
            }

            branch_div.attr('id', 'regression-div-' + i);
            branch_div.hide();
            main_div.append(branch_div);

            create_data_table(display_table, ignored_table, ignored_conf_sample_div,
                              data, params, branch, all_ignored_keys);
            branch_div.append(display_table);
            ignored_table.hide();
            ignored_conf_sample_div.hide();

            branch_div.append(ignored_table);
            branch_div.append(ignored_conf_sample_div);

            update_ignore_conf_sample(data, ignored_conf_sample_div, branch);

            branch_div.append(ignored_button);
            ignored_button.on('click', function(evt) {
                ignored_button.hide();
                $(".ignored").show();
            });
        });

        if (branches && branches.length > 1) {
            current_title = "Regressions in " + branches[0] + " branch";
        }
        $("#title").text(current_title);
        main_div.find("#regression-div-0").show();
        main_div.show();

        if (local_storage_available) {
            /* Clear out local storage space */
            var keys = Object.keys(localStorage);
            $.each(keys, function(i, key) {
                if (key.slice(0, ignore_key_prefix.length) == ignore_key_prefix &&
                        !all_ignored_keys[key]) {
                    delete localStorage[key];
                }
            });
        }

        return main_div;
    }

    function create_data_table(display_table, ignored_table, ignored_conf_sample_div,
                               data, params, branch, all_ignored_keys) {
        var table_head = $('<thead><tr>' +
                           '<th data-sort="string">Benchmark</th>' +
                           '<th data-sort="string">Date</th>' +
                           '<th data-sort="string">Commit(s)</th>' +
                           '<th data-sort="factor">Factor</th>' +
                           '<th data-sort="value">Best</th>' +
                           '<th data-sort="value">Current</th>' +
                           '<th></th>' +
                           '</tr></thead>');

        display_table.append(table_head);
        ignored_table.append(table_head.clone());

        var table_body = $('<tbody/>');
        var ignored_table_body = $('<tbody/>');

        var regressions = data['regressions'];

        $.each(regressions, function (i, item) {
            var benchmark_name = item[0];
            var graph_url = item[1];
            var param_dict = item[2];
            var parameter_idx = item[3];
            var regression = item[4];

            if (regression === null) {
                return;
            }

            if (branch !== null && param_dict['branch'] != branch) {
                return;
            }
            var revisions = regression[0];
            var new_value = regression[1];
            var old_value = regression[2];

            var factor = new_value / old_value;
            var last_revision = revisions[revisions.length - 1];
            var date_fmt = new Date($.asv.master_json.revision_to_date[last_revision[1]]);

            var row = $('<tr/>');

            var benchmark_basename = benchmark_name.replace(/\(.*/, '');
            var url_params = {};

            $.each(param_dict, function (key, value) {
                url_params[key] = [value];
            });

            url_params.commits = [];
            $.each(revisions, function(i, revs) {
                var commit_b = $.asv.get_commit_hash(revs[1]);
                if (revs[0] !== null) {
                    var commit_a = $.asv.get_commit_hash(revs[0]);
                    url_params.commits.push(commit_a + '-' + commit_b);
                }
                else {
                    url_params.commits.push(commit_b);
                }
            });

            if (parameter_idx !== null) {
                url_params.idx = [parameter_idx];
            }
            var benchmark_url = $.asv.format_hash_string({
                location: [benchmark_basename],
                params: url_params
            });

            if ($.asv.master_json.benchmarks[benchmark_basename].unit == "seconds") {
                new_value = $.asv.pretty_second(new_value);
                old_value = $.asv.pretty_second(old_value);
            }
            else {
                new_value = new_value.toPrecision(3);
                old_value = old_value.toPrecision(3);
            }

            var benchmark_link = $('<a/>').attr('href', benchmark_url).text(benchmark_name);
            row.append($('<td/>').append(benchmark_link));
            row.append($('<td class="date"/>').text($.asv.format_date_yyyymmdd_hhmm(date_fmt)));

            var commit_td = $('<td/>');
            $.each(revisions, function(i, revs) {
                var commit_a = $.asv.get_commit_hash(revs[0]);
                var commit_b = $.asv.get_commit_hash(revs[1]);
                if (i > 0) {
                    commit_td.append($('<span>, </span>'));
                }
                if (commit_a) {
                    if ($.asv.master_json.show_commit_url.match(/.*\/\/github.com\//)) {
                        var commit_url = ($.asv.master_json.show_commit_url + '../compare/'
                                          + commit_a + '...' + commit_b);
                        commit_td.append(
                            $('<a/>').attr('href', commit_url).text(commit_a + '..' + commit_b));
                    }
                    else {
                        commit_td.append($('<span/>').text(commit_a + '..' + commit_b));
                    }
                }
                else {
                    var commit_url = $.asv.master_json.show_commit_url + commit_b;
                    commit_td.append(
                        $('<a/>').attr('href', commit_url).text(commit_b));
                }
            });
            row.append(commit_td);

            row.append($('<td/>').text(factor.toFixed(2) + 'x'));
            row.append($('<td/>').text(old_value));
            row.append($('<td/>').text(new_value));

            /* html5 local storage has limited size, so store hashes
               rather than potentially long strings */
            var ignore_key = get_ignore_key(item);
            all_ignored_keys[ignore_key] = 1;

            var is_ignored = is_key_ignored(ignore_key);
            var ignore_button = $('<button class="btn btn-small"/>');

            row.attr('id', ignore_key);

            ignore_button.on('click', function(evt) {
                if (is_key_ignored(ignore_key)) {
                    set_key_ignore_status(ignore_key, false);
                    var item = ignored_table_body.find('#' + ignore_key).detach();
                    ignore_button.text('Ignore');
                    table_body.append(item);
                }
                else {
                    set_key_ignore_status(ignore_key, true);
                    var item = table_body.find('#' + ignore_key).detach();
                    ignore_button.text('Unignore');
                    ignored_table_body.append(item);
                }
                update_ignore_conf_sample(data, ignored_conf_sample_div, branch);
            });

            row.append($('<td/>').append(ignore_button));

            if (!is_ignored) {
                ignore_button.text('Ignore');
                table_body.append(row);
            }
            else {
                ignore_button.text('Unignore');
                ignored_table_body.append(row);
            }

            /* Show a graph as a popup */
            $.asv.ui.hover_graph(benchmark_link, graph_url, benchmark_basename, parameter_idx, revisions);
        });

        display_table.append(table_body);
        ignored_table.append(ignored_table_body);

        setup_sort(params, display_table);
        setup_sort(params, ignored_table);
    }

    function get_ignore_key(item) {
        var benchmark_name = item[0];
        var regression = item[4];

        if (regression === null) {
            return null;
        }
        var revisions = regression[0];
        var ignore_payload = benchmark_name;

        $.each(revisions, function (i, revs) {
            if (revs[0] === null) {
                ignore_payload = ignore_payload + ',';
            }
            else {
                ignore_payload = (ignore_payload + ','
                                  + $.asv.master_json.revision_to_hash[revs[0]]);
            }
            ignore_payload = (ignore_payload + ','
                              + $.asv.master_json.revision_to_hash[revs[1]]);
        });

        return ignore_key_prefix + $.md5(ignore_payload);
    }

    function is_key_ignored(ignore_key) {
        if (local_storage_available) {
            return (ignore_key in localStorage) || (ignore_key in ignored_regressions);
        }
        else {
            return (ignore_key in ignored_regressions);
        }
    }

    function set_key_ignore_status(ignore_key, is_ignored) {
        if (is_ignored) {
            if (local_storage_available) {
                try {
                    localStorage[ignore_key] = 1;
                } catch (err) {
                    /* Out of quota -- we're just going to ignore that */
                }
            }
            ignored_regressions[ignore_key] = 1;
        }
        else {
            if (local_storage_available) {
                delete localStorage[ignore_key];
            }
            delete ignored_regressions[ignore_key];
        }
    }

    function update_ignore_conf_sample(data, ignored_conf_sample_div, branch) {
        var regressions = data['regressions'];
        var entries = {};
        var branch_suffix = "";

        if (branch) {
            branch_suffix = "@" + branch;
        }

        $.each(regressions, function (i, item) {
            var param_dict = item[2];
            if (branch !== null && param_dict['branch'] != branch) {
                return;
            }

            var ignore_key = get_ignore_key(item);
            if (is_key_ignored(ignore_key)) {
                var benchmark_name = item[0];
                var regression = item[4];
                if (regression === null) {
                    return;
                }
                var benchmark_name_re = (benchmark_name + branch_suffix).replace(/[.?*+^$[\]\\(){}|-]/g, "\\\\$&");
                var revisions = regression[0];
                var last_commit = $.asv.get_commit_hash(revisions[revisions.length-1][1]);
                var entry = "        \"^" + benchmark_name_re + "$\": \"" + last_commit + "\",\n";
                entries[entry] = 1;
            }
        });

        entries = Object.keys(entries);
        entries.sort();

        var text = "// asv.conf.json excerpt for ignoring the above permanently\n\n";
        text += "    \"regressions_first_commits\": {\n";
        $.each(entries, function (i, entry) {
            text += entry;
        });
        text += "    }";

        var pre = $('<pre/>');
        pre.text(text);
        ignored_conf_sample_div.empty();
        ignored_conf_sample_div.append(pre);
    }

    function setup_sort(params, table) {
        table.stupidtable({
            'value': function(a, b) {
                function key(s) {
                    for (var k = 0; k < $.asv.time_units.length; ++k) {
                        var entry = $.asv.time_units[k];
                        var m = s.match('^([0-9.]+)'+entry[0]+'$');
                        if (m) {
                            return parseFloat(m[1]) * entry[2] * 1e-30;
                        }
                    }
                    return 0;
                }
                return key(a) - key(b)
            },
            'factor': function(a, b) {
                return parseFloat(a.replace(/x/, '')) - parseFloat(b.replace(/x/, ''));
            }
        });

        table.bind('aftertablesort', function (event, data) {
            var info = $.asv.parse_hash_string(window.location.hash);
            info.params['sort'] = [data.column];
            info.params['dir'] = [data.direction];
            window.location.hash = $.asv.format_hash_string(info);

            /* Update appearance */
            table.find('thead th').removeClass('asc');
            table.find('thead th').removeClass('desc');
            var th_to_sort = table.find("thead th").eq(parseInt(data.column));
            if (th_to_sort) {
                th_to_sort.addClass(data.direction);
            }
        });

        if (params.sort && params.dir) {
            var th_to_sort = table.find("thead th").eq(parseInt(params.sort[0]));
            th_to_sort.stupidsort(params.dir[0]);
        }
        else {
            var th_to_sort = table.find("thead th").eq(3);
            th_to_sort.stupidsort("desc");
        }
    }

    /*
      Setup display hooks
    */
    $.asv.register_page('regressions', function(params) {
        $('#regressions-display').show()
        load_data(params);
    });
});
