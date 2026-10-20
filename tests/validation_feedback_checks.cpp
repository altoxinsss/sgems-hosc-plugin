#include "../shared/ValidationFeedback.h"

#include <qcombobox.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <iostream>

namespace hv = hosc_validation;

static int failures = 0;

static void check(bool condition, const char* description)
{
    if(!condition) {
        ++failures;
        std::cerr << "FAIL: " << description << std::endl;
    }
}

static void test_error_messages()
{
    Error_messages_handler_xml errors;
    errors.report(true, "template_sampling",
        "legacy_search cannot sample a 3D grid; choose exhaustive_integer_offsets.");
    errors.report(true, "integer_direction_vectors",
        "Provide three vectors: dx,dy,dz separated by semicolons.");
    std::vector<hv::Issue> issues = hv::parse_issues(errors.errors());
    check(issues.size() == 2, "all field errors returned by the real S-GeMS handler are parsed");
    if(issues.size() == 2) {
        check(issues[0].field == "template_sampling", "first field key is preserved");
        check(issues[0].message.find("legacy_search") != std::string::npos,
            "first field explanation is preserved");
        check(issues[1].field == "integer_direction_vectors", "second field key is preserved");
    }
    issues = hv::parse_issues(
        "\n<field diagnostic=\"Use &lt;3 &amp; &quot;quoted&quot; values\" />\n"
        "<other diagnostic=\"Old unescaped \"quoted\" text\" />\n"
        "A general input error\n");
    check(issues.size() == 3, "blank lines are ignored while generic errors remain visible");
    if(issues.size() == 3) {
        check(issues[0].message == "Use <3 & \"quoted\" values", "escaped diagnostic text is decoded");
        check(issues[1].message == "Old unescaped \"quoted\" text", "legacy unescaped quotes are readable");
        check(issues[2].message == "A general input error", "non-field diagnostic text is preserved");
    }
    check(hv::parse_issues("").empty(), "empty errors produce an empty result");
}

static void test_hidden_form()
{
    QWidget form(0, "input_parameters");
    QTextEdit* summary = new QTextEdit(&form, "hosc_template_validation_summary");
    summary->setReadOnly(true);
    QTabWidget* tabs = new QTabWidget(&form, "parameter_tabs");
    QWidget* general = new QWidget(tabs, "general_page");
    QWidget* exact = new QWidget(tabs, "exact_template_page");
    tabs->addTab(general, "General");
    tabs->addTab(exact, "Exact Template");
    QScrollView* scroll = new QScrollView(exact, "exact_scroll");
    QWidget* controls = new QWidget(scroll->viewport(), "exact_controls");
    scroll->addChild(controls);
    QLabel* label = new QLabel("Sampling mode", controls, "sampling_label");
    QComboBox* sampling = new QComboBox(false, controls, "template_sampling");
    sampling->insertItem("legacy_search");
    sampling->insertItem("exhaustive_integer_offsets");
    QLineEdit* directions = new QLineEdit(controls, "integer_direction_vectors");
    directions->setText("1,0;0,1");
    QToolTip::add(sampling, "Choose how templates are sampled.");
    QWhatsThis::add(directions, "One integer vector per template direction.");
    tabs->setCurrentPage(0);
    const QPalette samplingPalette = sampling->palette();
    const QPalette directionsPalette = directions->palette();
    const QPalette labelPalette = label->palette();

    // Command-file execution uses the ordinary handler. It must not alter
    // an unrelated parameter form, even if a matching marker exists.
    summary->setText("Batch execution must leave this form untouched.");
    {
        Error_messages_handler batchErrors;
        hv::FeedbackGuard guard(&batchErrors, "hosc_template_validation_summary");
        batchErrors.report("template_sampling", "A batch-only failure.");
    }
    check(summary->text() == "Batch execution must leave this form untouched.",
        "ordinary batch error handler does not alter the GUI summary");
    check(sampling->palette() == samplingPalette, "ordinary batch handler does not highlight GUI controls");

    hv::reset_form(&form, summary);
    check(summary->text().contains("Ready."), "initial summary tells the user how to validate");
    Error_messages_handler_xml errors;
    errors.report(true, "template_sampling", "3D input requires exhaustive_integer_offsets.");
    errors.report(true, "integer_direction_vectors", "Expected three 3D integer vectors.");
    const std::vector<hv::Issue> issues = hv::parse_issues(errors.errors());
    hv::present(&form, summary, issues);
    check(summary->text().contains("2 input problem(s)"), "summary reports both input problems");
    check(summary->text().contains("Sampling mode") && summary->text().contains("Integer direction vectors"),
        "summary names the controls users see");
    check(summary->text().contains("3D input requires") && summary->text().contains("Expected three"),
        "summary includes concrete corrections");
    check(summary->textFormat() == Qt::PlainText, "diagnostics cannot be interpreted as rich text markup");
    check(tabs->currentPageIndex() == 1, "first failing control's tab is selected");
    const QColor errorColor(255, 205, 205);
    check(sampling->palette().active().button() == errorColor,
        "non-editable combo uses red Button role");
    check(sampling->palette().inactive().button() == errorColor,
        "combo remains red when the host displays its dialog");
    check(directions->palette().active().base() == errorColor,
        "editable text control uses red Base role");
    check(label->palette() == labelPalette, "unrelated label palette is unchanged");
    check(QToolTip::textFor(sampling).contains("3D input requires"), "field tooltip explains the failure");
    check(!form.isVisible(), "validation does not display the hidden test form");

    // A repeated failure must not replace the original palette with red.
    hv::present(&form, summary, issues);
    hv::reset_form(&form, summary);
    check(sampling->palette() == samplingPalette, "reset restores the original combo palette");
    check(directions->palette() == directionsPalette, "reset restores the original text palette");
    check(QToolTip::textFor(sampling) == "Choose how templates are sampled.", "reset restores original tooltip");
    check(QWhatsThis::textFor(directions) == "One integer vector per template direction.",
        "reset restores original help text");
    check(!sampling->child("hosc_saved_validation_palette", "QObject", false),
        "reset removes saved palette state");
    const std::vector<hv::Issue> none;
    hv::present(&form, summary, none);
    check(summary->text().contains("Parameter validation passed."), "success replaces stale failure text");
    check(!summary->text().contains("input problem(s)"), "success summary has no stale error count");

    hv::Issue unknown;
    unknown.field = "not_a_widget";
    unknown.message = "A general check failed.";
    std::vector<hv::Issue> unknownIssues(1, unknown);
    hv::present(&form, summary, unknownIssues);
    check(summary->text().contains("A general check failed."), "unmapped errors are still shown in summary");
}

int main(int argc, char** argv)
{
    QApplication application(argc, argv);
    test_error_messages();
    test_hidden_form();
    if(failures) {
        std::cerr << failures << " validation-feedback check(s) failed." << std::endl;
        return 1;
    }
    std::cout << "Validation-feedback checks passed: error parsing, full summary, tab selection, "
        "combo/text highlighting, tooltip preservation, repeated-failure reset and success state."
        << std::endl;
    return 0;
}
