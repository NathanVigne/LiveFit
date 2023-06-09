// Copyright (C) 2021 The Qt Company Ltd.
// SPDX-License-Identifier: LicenseRef-Qt-Commercial OR GPL-3.0-only

#ifndef MAINWIDGET_H
#define MAINWIDGET_H

#include <QSlider>
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QScatterSeries>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QWidget>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>

#include <QLabel>
#include "worker.h"

QT_USE_NAMESPACE

class MainWidget : public QWidget
{
    Q_OBJECT
public:
    explicit MainWidget(QWidget *parent = nullptr);
    void createUi();

signals:
    void seriesChange(double *datas);
    void drawFit();

public Q_SLOTS:
    void startFitting();
    void getSliderValue();
    void copyParam(double *fit_params);
    void drawingFit();

protected:
    void updataSeries();
    void createSeries();

private:
    QChart *m_chart;
    QChartView *m_chartView;
    QScatterSeries *m_serieData;
    QLineSeries *m_serieFit;

    worker *work;
    double *datas;
    double param[5];
    double fitParam[4];
    size_t N;
    std::mutex mutex;

    QPushButton *m_startFit;
    QSlider *A_slider;
    QLabel *A_label;
    QSlider *x_slider;
    QLabel *x_label;
    QSlider *w_slider;
    QLabel *w_label;
    QSlider *b_slider;
    QLabel *b_label;
    QSlider *noise_slider;
    QLabel *noise_label;
};

#endif // MAINWIDGET_H
