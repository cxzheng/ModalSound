/******************************************************************************
 *  File: ParamsDialog.h
 *
 *  This file is part of isostuffer
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
#ifndef UI_PARAMS_DIALOG_H
#   define UI_PARAMS_DIALOG_H

#include <QDoubleValidator>

#if QT_VERSION >= 0x040000
#include "ui_Params.h"
#else
#error Only Qt with the version higher than 4.0 is supported right now
#endif

class IsoStufferFrame;

class ParamsDialog : public QDialog, private Ui_ParamsDialog
{
    friend class IsoStufferFrame;

    Q_OBJECT

    public slots:
        void res_changed(int res)
        {   
            levelSpin->setMaximum(res+1); 
            marginSpin->setMaximum(1 << res);
        }

    public:
        ParamsDialog():validator_(this)
        {
            setupUi(this);
            alphaLongBox->setValidator(&validator_);
            alphaShortBox->setValidator(&validator_);

            QObject::connect(resSpin, SIGNAL(valueChanged(int)), this, SLOT(res_changed(int)));
        }

        double alpha_short() const
        {   return alphaShortBox->text().toDouble(); }

        double alpha_long() const
        {   return alphaLongBox->text().toDouble(); }

        int resolution() const
        {   return (1 << resSpin->value()); }

        int num_levels() const
        {   return levelSpin->value(); }

        int margin() const
        {   return marginSpin->value(); }

    private:
        QDoubleValidator    validator_;
};

#endif
