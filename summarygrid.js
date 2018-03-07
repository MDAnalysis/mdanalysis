'use strict';

$(document).ready(function() {
    var summary_loaded = false;

    /* Callback a function when an element comes in view */
    function callback_in_view(element, func) {
        function handler(evt) {
            var visible = (
                $('#summarygrid-display').css('display') != 'none' &&
                (element.offset().top <= $(window).height() + $(window).scrollTop()) &&
                    (element.offset().top + element.height() >= $(window).scrollTop()));
            if (visible) {
                func();
                $(window).unbind('scroll', handler);
            }
        }
        $(window).on('scroll', handler);
    }

    function make_summary() {
        var summary_display = $('#summarygrid-display');
        var master_json = $.asv.master_json;
        var summary_container = $('<div/>');

        if (summary_loaded) {
            return;
        }

        $.each(master_json.benchmarks, function(bm_name, bm) {
            var container = $(
                '<a class="btn benchmark-container" href="#' + bm_name +
                '"/>');
            var plot_div = $(
                '<div id="summarygrid-' + bm_name + '" class="benchmark-plot"/>');
            var name = $('<div class="benchmark-text">' + bm_name + '</div>');
            name.tooltip({
                title: bm_name,
                html: true,
                placement: 'top',
                container: 'body',
                animation: false
            });

            plot_div.tooltip({
                title: bm.code,
                html: true,
                placement: 'bottom',
                container: 'body',
                animation: false
            });

            container.append(name);
            container.append(plot_div);
            summary_container.append(container);

            callback_in_view(plot_div, function() {
                $.asv.load_graph_data(
                    'graphs/summary/' + bm_name + '.json'
                ).done(function(data) {
                    var options = {
                        colors: $.asv.colors,
                        series: {
                            lines: {
                                show: true,
                                lineWidth: 2
                            },
                            shadowSize: 0
                        },
                        grid: {
                            borderWidth: 1,
                            margin: 0,
                            labelMargin: 0,
                            axisMargin: 0,
                            minBorderMargin: 0
                        },
                        xaxis: {
                            ticks: [],
                        },
                        yaxis: {
                            ticks: [],
                            min: 0
                        },
                        legend: {
                            show: false
                        }
                    };

                    var plot = $.plot(
                        plot_div, [{data: data}], options);
                }).fail(function() {
                    // TODO: Handle failure
                });
            });
        });

        summary_display.append(summary_container);
        $(window).scroll();

        summary_loaded = true;
    }

    $.asv.register_page('', function(params) {
        $('#summarygrid-display').show();
        $("#title").text("All benchmarks");
        $('.tooltip').remove();
        make_summary();
    });
});
