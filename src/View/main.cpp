#include "view.h"
#include "widget.h"
//#include "model.h"
//#include "controller.h"

#include <QApplication>

int main(int argc, char *argv[]) {
  QApplication a(argc, argv);
  s21::Model m;
  s21::GraphModel gm;
  s21::Controller c(&m, &gm);
  View v(&c);
  v.show();

  return a.exec();
}
