#ifndef UI_MODES_VIS_PARAMS_DIAG_INC
#   define UI_MODES_VIS_PARAMS_DIAG_INC

#include <QDoubleValidator>                                    
                                                               
#if QT_VERSION >= 0x040000                                     
#include "ui_ModesVisParams.h"                                    
#else                                                          
#error Only Qt with the version higher than 4.0 is supported right now
#endif  

class ModesVisParamsDialog : public QDialog, private Ui_ModesVisParamsDialog
{
    Q_OBJECT

    signals:
        void scaleChanged(double);

    public slots:
        void scale_changed()
        {   emit scaleChanged(boxModeScale_->text().toDouble()); }

    public:
        ModesVisParamsDialog(QObject* parent):validator_(this)
        {
            setupUi(this);
            boxModeScale_->setValidator(&validator_);
            QObject::connect(boxModeId_, SIGNAL(valueChanged(int)), parent, SLOT(set_mode_sel(int)));
            QObject::connect(boxModeScale_, SIGNAL(editingFinished()), this, SLOT(scale_changed()));
        }

        void set_mode_num(int n)
        {
            boxModeId_->setMinimum(0);
            boxModeId_->setMaximum(n-1);
        }

    private:
        QDoubleValidator    validator_;
};
                                                               
#endif
