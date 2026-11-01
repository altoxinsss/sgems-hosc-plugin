#ifndef HOSC_VALIDATION_FEEDBACK_H
#define HOSC_VALIDATION_FEEDBACK_H

#include <GsTLAppli/utils/error_messages_handler.h>
#include <qapplication.h>
#include <qobjectlist.h>
#include <qwidgetlist.h>
#include <qtextedit.h>
#include <qtabwidget.h>
#include <qscrollview.h>
#include <qtooltip.h>
#include <qwhatsthis.h>
#include <sstream>
#include <string>
#include <vector>

// Qt3-only presentation adapter. Numerical validation remains in the plugins.
// No dialog or event loop is started here: scripted runs must never wait for UI.
namespace hosc_validation {

struct Issue {
    std::string field;
    std::string message;
};

inline void replace_entity(std::string& text, const char* entity, const char* value)
{
    const std::string key(entity);
    size_t offset = 0;
    while((offset = text.find(key, offset)) != std::string::npos) {
        text.replace(offset, key.size(), value);
        offset += std::string(value).size();
    }
}

inline std::string decode_message(std::string text)
{
    replace_entity(text, "&quot;", "\"");
    replace_entity(text, "&apos;", "'");
    replace_entity(text, "&lt;", "<");
    replace_entity(text, "&gt;", ">");
    replace_entity(text, "&#10;", "\n");
    replace_entity(text, "&#13;", "\r");
    replace_entity(text, "&amp;", "&");
    return text;
}

inline std::vector<Issue> parse_issues(const std::string& text)
{
    // Error_messages_handler emits one <field diagnostic="..." /> per line.
    // Use the last quote so an old, unescaped diagnostic is still readable.
    std::vector<Issue> issues;
    std::istringstream input(text);
    std::string line;
    while(std::getline(input, line)) {
        if(line.find_first_not_of(" \t\r") == std::string::npos) continue;
        Issue issue;
        const size_t start = line.find('<');
        const size_t marker = line.find("diagnostic=\"");
        const size_t last_quote = line.rfind('"');
        if(start != std::string::npos && marker != std::string::npos &&
           marker > start && last_quote >= marker + 12) {
            const size_t end = line.find_first_of(" \t>", start + 1);
            issue.field = line.substr(start + 1, end - start - 1);
            issue.message = decode_message(line.substr(marker + 12, last_quote - marker - 12));
        }
        else {
            issue.field = "general";
            issue.message = line;
        }
        issues.push_back(issue);
    }
    return issues;
}

class SavedPalette : public QObject {
public:
    SavedPalette(QWidget* widget)
        : QObject(widget, "hosc_saved_validation_palette"), palette(widget->palette()),
          tooltip(QToolTip::textFor(widget)), whatsThis(QWhatsThis::textFor(widget)) {}
    QPalette palette;
    QString tooltip;
    QString whatsThis;
};

inline void reset_form(QWidget* form, QTextEdit* summary)
{
    QObjectList* states = form->queryList("QObject", "hosc_saved_validation_palette", false, true);
    for(QObject* state = states->first(); state; state = states->next()) {
        SavedPalette* saved = static_cast<SavedPalette*>(state);
        QWidget* widget = static_cast<QWidget*>(saved->parent());
        widget->setPalette(saved->palette);
        QToolTip::remove(widget);
        QWhatsThis::remove(widget);
        if(!saved->tooltip.isEmpty()) QToolTip::add(widget, saved->tooltip);
        if(!saved->whatsThis.isEmpty()) QWhatsThis::add(widget, saved->whatsThis);
    }
    // queryList returns a separate list; deleting its entries is safe here.
    states->setAutoDelete(true);
    delete states;
    summary->setTextFormat(Qt::PlainText);
    summary->setText("Ready. Click Run to validate the parameters.");
}

inline QString field_title(QWidget* field, const std::string& name)
{
    const char* label = name.c_str();
    if(name == "Source_Grid") label = "Source Grid";
    else if(name == "Source_Property") label = "Source Property";
    else if(name == "props") label = "Target Property";
    else if(name == "Target_Button") label = "Calculate Target";
    else if(name == "Algorithm_Mode" || name == "Algorithm_Button") label = "Calculate Algorithm";
    else if(name == "template_sampling") label = "Sampling mode";
    else if(name == "integer_direction_vectors") label = "Integer direction vectors";
    else if(name == "Weighted_Enable") label = "Parameter file / template weighting";
    else if(name == "normalize_weights") label = "Parameter file / weight normalization";
    else if(name == "zero_first_order_moment") label = "Parameter file / legacy centering";
    else if(name == "invert_window_output") label = "Parameter file / primary-statistic inversion";
    else if(name == "weight_center_mode") label = "Parameter file / coordinate weight center";
    else if(name == "NOD") label = "Number of occurrences";
    else if(name == "NLAG") label = "Number of lags";
    else if(name == "LAG") label = "Lag separation";
    else if(name == "moment_order") label = "Moment/Cumulant Order";
    else if(name == "windowl") label = "Window Size X";
    else if(name == "windoww") label = "Window Size Y";
    else if(name == "windowh") label = "Window Size Z";
    QString title = QString::fromLocal8Bit(label);
    for(QWidget* ancestor = field; ancestor; ancestor = ancestor->parentWidget()) {
        QWidget* parent = ancestor->parentWidget();
        if(parent && parent->inherits("QTabWidget")) {
            QTabWidget* tabs = static_cast<QTabWidget*>(parent);
            // Pages normally live inside a QWidgetStack, so inspect ancestors
            // against each page rather than requiring a direct tab parent.
            for(int i = 0; i < tabs->count(); ++i) {
                QWidget* page = tabs->page(i);
                for(QWidget* item = field; item; item = item->parentWidget()) {
                    if(item == page) return tabs->tabLabel(page) + " / " + title;
                }
            }
        }
    }
    return title;
}

inline void reveal_field(QWidget* form, QWidget* field)
{
    QObjectList* tab_widgets = form->queryList("QTabWidget");
    for(QObject* object = tab_widgets->first(); object; object = tab_widgets->next()) {
        QTabWidget* tabs = static_cast<QTabWidget*>(object);
        for(int i = 0; i < tabs->count(); ++i) {
            QWidget* page = tabs->page(i);
            for(QWidget* ancestor = field; ancestor; ancestor = ancestor->parentWidget()) {
                if(ancestor == page) { tabs->showPage(page); break; }
            }
        }
    }
    delete tab_widgets;
    // Retain the summary at the top; scroll only nested tab scroll views.
    for(QWidget* ancestor = field->parentWidget(); ancestor && ancestor != form;
        ancestor = ancestor->parentWidget()) {
        if(ancestor->inherits("QScrollView")) {
            QScrollView* scroll = static_cast<QScrollView*>(ancestor);
            const QPoint p = scroll->viewportToContents(field->mapTo(scroll->viewport(), QPoint(0, 0)));
            scroll->ensureVisible(p.x(), p.y(), 20, 30);
        }
    }
    field->setFocus();
}

inline void present(QWidget* form, QTextEdit* summary, const std::vector<Issue>& issues)
{
    if(issues.empty()) {
        summary->setText("Parameter validation passed.");
        return;
    }
    std::ostringstream heading;
    heading << issues.size() << " input problem(s). Correct the fields below, then click Run.\n";
    QString text = QString::fromLocal8Bit(heading.str().c_str());
    QWidget* first = NULL;
    for(size_t i = 0; i < issues.size(); ++i) {
        QWidget* field = static_cast<QWidget*>(form->child(issues[i].field.c_str(), "QWidget", true));
        text += "\n" + QString::number((int)i + 1) + ". " + field_title(field, issues[i].field) + "\n";
        text += QString::fromLocal8Bit(issues[i].message.c_str()) + "\n";
        if(!field) continue;
        if(!first) first = field;
        if(!field->child("hosc_saved_validation_palette", "QObject", false)) new SavedPalette(field);
        QPalette palette = field->palette();
        // Non-editable QComboBox uses Button, while text inputs use Base.
        // Set all groups so the highlight survives the host's modal popup.
        palette.setColor(QColorGroup::Base, QColor(255, 205, 205));
        palette.setColor(QColorGroup::Button, QColor(255, 205, 205));
        palette.setColor(QColorGroup::Text, QColor(0, 0, 0));
        palette.setColor(QColorGroup::ButtonText, QColor(0, 0, 0));
        field->setPalette(palette);
        QToolTip::add(field, QString::fromLocal8Bit(issues[i].message.c_str()));
        QWhatsThis::add(field, QString::fromLocal8Bit(issues[i].message.c_str()));
    }
    summary->setTextFormat(Qt::PlainText);
    summary->setText(text);
    summary->setCursorPosition(0, 0);
    if(first) reveal_field(form, first);
    // The host wraps the entire form in a scroll view. A previous visit to
    // the bottom of the form must not leave the new summary off-screen.
    for(QWidget* ancestor = form->parentWidget(); ancestor; ancestor = ancestor->parentWidget()) {
        if(ancestor->inherits("QScrollView")) {
            QScrollView* scroll = static_cast<QScrollView*>(ancestor);
            const QPoint p = scroll->viewportToContents(summary->mapTo(scroll->viewport(), QPoint(0, 0)));
            scroll->ensureVisible(p.x(), p.y(), 10, summary->height());
            break;
        }
    }
}

class FeedbackGuard {
public:
    FeedbackGuard(Error_messages_handler* errors, const char* summary_name)
        : errors_(errors), form_(NULL), summary_(NULL)
    {
        // Command files also create QApplication. The handler type identifies
        // the GUI Run action; never search or modify forms during batch runs.
        if(!qApp || !dynamic_cast<Error_messages_handler_xml*>(errors)) return;
        QWidgetList* widgets = QApplication::allWidgets();
        for(QWidget* widget = widgets->first(); widget; widget = widgets->next()) {
            if(std::string(widget->name()) == summary_name &&
               widget->inherits("QTextEdit") && widget->isVisible()) {
                summary_ = static_cast<QTextEdit*>(widget);
                form_ = widget->parentWidget();
                break;
            }
        }
        delete widgets;
        if(form_ && summary_) reset_form(form_, summary_);
    }
    ~FeedbackGuard()
    {
        // Validation feedback must not introduce a new exception boundary.
        try {
            if(form_ && summary_ && errors_) present(form_, summary_, parse_issues(errors_->errors()));
        }
        catch(...) {}
    }
private:
    FeedbackGuard(const FeedbackGuard&);
    FeedbackGuard& operator=(const FeedbackGuard&);
    Error_messages_handler* errors_;
    QWidget* form_;
    QTextEdit* summary_;
};

} // namespace hosc_validation
#endif
