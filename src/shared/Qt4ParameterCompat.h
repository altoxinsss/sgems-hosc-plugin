#ifndef HOSC_QT4_PARAMETER_COMPAT_H
#define HOSC_QT4_PARAMETER_COMPAT_H
#include <GsTLAppli/gui/QWidget_value_accessors/qcombobox_accessor.h>
#include <GsTLAppli/appli/manager_repository.h>
#include <QApplication>
#include <QVariant>
#include <QDomDocument>
#include <QTimerEvent>

namespace hosc_qt4 {
// The 2011 Beta accessor still invokes Qt3's removed internalActivate slot.
// Only HOSC forms use the replacement behavior; other forms retain the host path.
class ComboAccessor : public QComboBox_accessor {
    QComboBox* box_;
public:
    ComboAccessor():box_(0){}
    static Named_interface* create(std::string&){return new ComboAccessor;}
    virtual bool initialize(QWidget* widget=0){box_=qobject_cast<QComboBox*>(widget);return QComboBox_accessor::initialize(widget);}
    virtual bool set_value(const std::string& text){
        bool hosc=false;
        for(QWidget* p=box_;p;p=p->parentWidget())if(p->property("hoscQt4ParameterForm").toBool()){hosc=true;break;}
        if(!hosc)return QComboBox_accessor::set_value(text);
        QDomDocument doc;if(!doc.setContent(QString("<root>")+QString::fromStdString(text)+"</root>"))return false;
        QString wanted=doc.documentElement().firstChildElement().attribute("value");
        int index=box_->findText(wanted);
        if(index<0){if(!box_->isEditable())return false;box_->addItem(wanted);index=box_->count()-1;}
        box_->setCurrentIndex(index);
        QMetaObject::invokeMethod(box_,"activated",Qt::DirectConnection,Q_ARG(int,index));
        QMetaObject::invokeMethod(box_,"activated",Qt::DirectConnection,Q_ARG(QString,wanted));
        return true;
    }
};
class Installer : public QObject {
    int timer_;
public:
    explicit Installer(QObject* parent):QObject(parent),timer_(startTimer(0)){}
protected:
    virtual void timerEvent(QTimerEvent*){
        SmartPtr<Named_interface> ni=Root::instance()->interface(xmlGeostatParamUtils_manager);
        Manager* manager=dynamic_cast<Manager*>(ni.raw_ptr());if(!manager)return;
        manager->factory("HOSC_Qt4ComboAccessor",ComboAccessor::create);
        manager->delete_interface("/QComboBox");
        SmartPtr<Named_interface> added=manager->new_interface("HOSC_Qt4ComboAccessor","/QComboBox");
        if(added.raw_ptr()){killTimer(timer_);deleteLater();}
    }
};
inline void install(){
    if(!qApp||qApp->property("hoscQt4ParameterCompatInstalled").toBool())return;
    qApp->setProperty("hoscQt4ParameterCompatInstalled",true);
    new Installer(qApp);
}
}
#endif
