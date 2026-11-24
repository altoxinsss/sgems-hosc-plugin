#include <GsTLAppli/appli/library_appli_init.h>
#include <GsTLAppli/actions/library_actions_init.h>
#include <GsTLAppli/filters/library_filters_init.h>
#include <GsTLAppli/geostat/library_geostat_init.h>
#include <GsTLAppli/grid/library_grid_init.h>
#include <GsTLAppli/gui/library_gui_init.h>
#include <GsTLAppli/extra/gui/library_extragui_init.h>
#include <GsTLAppli/appli/project.h>
#include <GsTLAppli/appli/manager_repository.h>
#include <GsTLAppli/gui/utils/qwidget_value_collector.h>
#include <GsTLAppli/gui/appli/qt_algo_control_panel.h>
#include <GsTLAppli/utils/error_messages_handler.h>
#include <QtDesigner/QFormBuilder>
#include <QApplication>
#include <QFile>
#include <QTabWidget>
#include <QLineEdit>
#include <QComboBox>
#include <QVariant>
#include <QLibrary>
#include <QCheckBox>
#include <QPixmap>
#include <QFont>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include "../src/shared/ValidationFeedback.h"

static void require(bool good,const char* message) {
    if(!good)throw std::runtime_error(message);
}
template<class T> static T* field(QWidget* form,const char* name) {
    T* result=form->findChild<T*>(name);
    require(result!=0,name);return result;
}
int main(int argc,char**argv) {
  try {
    require(argc==4,"host, output, fixtures required");
    QApplication app(argc,argv);app.setFont(QFont("Arial",10));
    std::string host=argv[1],out=argv[2],fixtures=argv[3];
    Root::instance()->set_main_plugin_path(host+"/plugins/");
    std::cerr<<"BEGIN built-in host library initialization"<<std::endl;
    libGsTLAppli_actions_init();libGsTLAppli_appli_init();libGsTLAppli_filters_init();
    libGsTLAppli_geostat_init();libGsTLAppli_grid_init();libGsTLAppli_gui_init();libGsTLAppli_extragui_init();
    std::cerr<<"END built-in host library initialization; BEGIN HOSC checks"<<std::endl;
    QLibrary plugin((host+"/plugins/Geostat/HOSC2008.dll").c_str());
    typedef int (*PluginInit)();PluginInit init=(PluginInit)plugin.resolve("plugin_init");
    require(init&&init()==0,"register HOSC DLL");
    QLibrary window_plugin((host+"/plugins/Geostat/HOSCWindow.dll").c_str());
    PluginInit window_init=(PluginInit)window_plugin.resolve("plugin_init");
    require(window_init&&window_init()==0,"register window DLL");
    SmartPtr<Named_interface> template_probe=Root::instance()->new_interface("HOSC_TEMPLATE",geostatAlgo_manager+"/");
    require(template_probe.raw_ptr()!=0,"load HOSC DLL before UI");
    app.processEvents();
    GsTL_project project;Error_messages_handler_xml errors;
    require(project.execute("LoadObjectFromFile",fixtures+"/grid.sgems::sgems",&errors),"load grid");
    require(project.execute("LoadObjectFromFile",fixtures+"/points.sgems::sgems",&errors),"load points");
    QApplication::addLibraryPath((host+"/plugins").c_str());
    Algo_control_panel* panel=new Algo_control_panel(&project);
    panel->resize(660,1150);panel->move(-20000,-20000);panel->show();app.processEvents();
    QWidgets_values_collector_xml collector;
    const char* names[]={"HOSC_TEMPLATE","HOSC_WINDOW"};
    for(int a=0;a<2;++a) {
      std::string name=names[a];QFile file((host+"/plugins/Geostat/"+name+".ui").c_str());
      require(file.open(QIODevice::ReadOnly),"open form");
      QList<QTreeWidget*> trees=panel->findChildren<QTreeWidget*>();QTreeWidgetItem* selected=0;
      for(int i=0;i<trees.size();++i){QList<QTreeWidgetItem*> matches=trees[i]->findItems(name.c_str(),Qt::MatchExactly|Qt::MatchRecursive);if(!matches.isEmpty())selected=matches.first();}
      require(selected!=0,"algorithm registered in actual host selector");
      selected->treeWidget()->expandAll();selected->treeWidget()->setCurrentItem(selected);
      panel->prompt_for_parameters(Qt::LeftButton,selected,0);app.processEvents();
      QWidget* form=panel->findChild<QWidget*>("input_parameters");require(form!=0,"load form with actual SGeMS QUiLoader");
      QComboBox* grid=field<QComboBox>(form,"Source_Grid");
      QComboBox* prop=field<QComboBox>(form,"Source_Property");
      QComboBox* wp=field<QComboBox>(form,"weight_property");
      require(grid->inherits("GridSelectorBasic")&&prop->inherits("SinglePropertySelector")&&wp->inherits("SinglePropertySelector"),"custom selectors");
      grid->setCurrentIndex(grid->findText("grid"));
      QMetaObject::invokeMethod(grid,"activated",Qt::DirectConnection,Q_ARG(QString,QString("grid")));app.processEvents();
      require(prop->findText("phi")>=0&&wp->findText("phi")>=0,"source/property connection");
      prop->setCurrentIndex(prop->findText("phi"));wp->setCurrentIndex(wp->findText("phi"));
      QTabWidget* tabs=field<QTabWidget>(form,a?"tabWidget3":"parameter_tabs");
      require(tabs->count()==(a?3:6),"tab count");
      QComboBox* weights=field<QComboBox>(form,"weight_mode");require(weights->count()==5,"weight choices");
      for(int k=0;k<5;++k) {
        weights->setCurrentIndex(k);
        std::string saved=collector.widgets_values(form,name);
        weights->setCurrentIndex((k+1)%5);collector.set_widgets_values(saved,form);
        require(weights->currentIndex()==k,"weight round trip");
      }
      weights->setCurrentIndex(0);
      if(a) {
        const char* keys[]={"center_mode","neighborhood_shape","grid_bounds_mode"};
        for(int k=0;k<3;++k) {
          QComboBox* box=field<QComboBox>(form,keys[k]);require(box->count()==2,"window choices");
          box->setCurrentIndex(1);std::string saved=collector.widgets_values(form,name);
          box->setCurrentIndex(0);collector.set_widgets_values(saved,form);
          require(box->currentIndex()==1,"window round trip");box->setCurrentIndex(0);
        }
      } else {
        QCheckBox* on=field<QCheckBox>(form,"enable_weighted_moment");
        QWidget* settings=field<QWidget>(form,"template_weight_settings");
        on->setChecked(false);app.processEvents();require(settings->isHidden(),"weight switch off");
        on->setChecked(true);app.processEvents();require(!settings->isHidden(),"weight switch on");
      }
      std::ofstream xml((out+"/"+name+"_roundtrip.xml").c_str());xml<<collector.widgets_values(form,name);xml.close();
      form->resize(620,a?930:850);form->move(-20000,-20000);form->show();app.processEvents();
      for(int t=0;t<tabs->count();++t) {
        tabs->setCurrentIndex(t);app.processEvents();
        char suffix[32];sprintf(suffix,"_%d.png",t);
        require(QPixmap::grabWidget(form).save((out+"/"+name+suffix).c_str(),"PNG"),"save UI image");
      }
      if(a) {
        tabs->setCurrentIndex(1);QTabWidget* inner=field<QTabWidget>(form,"pointset_layout_tabs");
        for(int t=0;t<inner->count();++t) {
          inner->setCurrentIndex(t);app.processEvents();char suffix[32];sprintf(suffix,"_inner_%d.png",t);
          require(QPixmap::grabWidget(form).save((out+"/"+name+suffix).c_str(),"PNG"),"save inner tab");
        }
      }
      require(QPixmap::grabWidget(panel).save((out+"/"+name+"_host_panel.png").c_str(),"PNG"),"save actual host panel");
      std::cout<<name<<": actual SGeMS selector/QUiLoader, Qt4 custom widgets, signals, all weight modes and parameter round trips passed"<<std::endl;
      panel->load_parameters((out+"/"+name+"_run.par").c_str());app.processEvents();
      panel->run_geostat_algo();app.processEvents();
      std::string object=name+"_ui_smoke_"+name;
      require(project.execute("SaveGeostatGrid",object+"::"+out+"/"+name+"_ui_result.sgems::sgems",&errors),"actual panel Run result");
      std::cout<<name<<": actual host Load/Run/Save passed"<<std::endl;
      QTextEdit* summary=field<QTextEdit>(form,a?"hosc_window_validation_summary":"hosc_template_validation_summary");
      QWidget* badfield=field<QWidget>(form,a?"moment_order":"NOD");
      QPalette old=badfield->palette();QString oldtip=badfield->toolTip();
      hosc_validation::Issue issue;issue.field=a?"moment_order":"NOD";issue.message="Test validation message <not markup>.";
      std::vector<hosc_validation::Issue> issues(1,issue);
      hosc_validation::present(form,summary,issues);
      require(summary->toPlainText().contains("<not markup>"),"plain text diagnostics");
      require(badfield->palette().color(QPalette::Base)==QColor(255,205,205),"error field highlight");
      hosc_validation::reset_form(form,summary);
      require(badfield->palette()==old&&badfield->toolTip()==oldtip,"reset error highlight and tooltip");
      std::cout<<name<<": Qt4 validation feedback/reset passed"<<std::endl;
    }
    delete static_cast<QWidget*>(panel);
    return 0;
  } catch(const std::exception& e) {std::cerr<<e.what()<<std::endl;return 1;}
}
