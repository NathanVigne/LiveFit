// Copyright (C) 2021 The Qt Company Ltd.
// SPDX-License-Identifier: LicenseRef-Qt-Commercial OR GPL-3.0-only

#include "mainwidget.h"
#include <QLabel>
#include <QValueAxis>
#include <QtCharts/QBarSeries>
#include <QtCharts/QBarSet>
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLegend>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QFormLayout>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QPushButton>
#include <random>
QT_USE_NAMESPACE

MainWidget::MainWidget(QWidget *parent)
    : QWidget(parent)
{
    // Create chart view with the chart
    m_chart = new QChart();
    m_chartView = new QChartView(m_chart, this);

    m_chart->setTitle("Test fit Gaussian");

    createUi();

    connect(m_startFit, &QPushButton::clicked, this, &MainWidget::startFitting);
    connect(this, &MainWidget::drawFit, this, &MainWidget::drawingFit);

    createSeries();

    QValueAxis *axisX = new QValueAxis;
    axisX->setRange(0, 1000);
    axisX->setTickCount(10);
    m_chart->addAxis(axisX, Qt::AlignBottom);
    m_serieData->attachAxis(axisX);
    m_serieFit->attachAxis(axisX);

    QValueAxis *axisY = new QValueAxis;
    axisY->setRange(0, 1000);
    axisY->setTickCount(10);
    m_chart->addAxis(axisY, Qt::AlignLeft);
    m_serieData->attachAxis(axisY);
    m_serieFit->attachAxis(axisY);

    m_chartView->setRenderHint(QPainter::Antialiasing);
}

void MainWidget::createUi()
{
    // Create buttons for ui
    m_startFit = new QPushButton("Start Fit");

    // Creates slider for param
    A_slider = new QSlider(Qt::Horizontal);
    A_slider->setRange(0, 1000);
    A_slider->setValue(850);
    A_label = new QLabel("A");
    QHBoxLayout *A_layout = new QHBoxLayout();
    A_layout->addWidget(A_label);
    A_layout->addWidget(A_slider);
    connect(A_slider, &QSlider::valueChanged, this, &MainWidget::updataSeries);

    x_slider = new QSlider(Qt::Horizontal);
    x_slider->setRange(0, 1000);
    x_slider->setValue(500);
    x_label = new QLabel("x0");
    QHBoxLayout *x_layout = new QHBoxLayout();
    x_layout->addWidget(x_label);
    x_layout->addWidget(x_slider);
    connect(x_slider, &QSlider::valueChanged, this, &MainWidget::updataSeries);

    w_slider = new QSlider(Qt::Horizontal);
    w_slider->setRange(0, 1000);
    w_slider->setValue(200);
    w_label = new QLabel("w0");
    QHBoxLayout *w_layout = new QHBoxLayout();
    w_layout->addWidget(w_label);
    w_layout->addWidget(w_slider);
    connect(w_slider, &QSlider::valueChanged, this, &MainWidget::updataSeries);

    b_slider = new QSlider(Qt::Horizontal);
    b_slider->setRange(0, 250);
    b_slider->setValue(50);
    b_label = new QLabel("b");
    QHBoxLayout *b_layout = new QHBoxLayout();
    b_layout->addWidget(b_label);
    b_layout->addWidget(b_slider);
    connect(b_slider, &QSlider::valueChanged, this, &MainWidget::updataSeries);

    noise_slider = new QSlider(Qt::Horizontal);
    noise_slider->setRange(0, 100);
    noise_slider->setValue(10);
    noise_label = new QLabel("b");
    QHBoxLayout *noise_layout = new QHBoxLayout();
    noise_layout->addWidget(noise_label);
    noise_layout->addWidget(noise_slider);
    connect(noise_slider, &QSlider::valueChanged, this, &MainWidget::updataSeries);

    // Create layout for grid and detached legend
    QVBoxLayout *control_Layout = new QVBoxLayout();
    control_Layout->addLayout(A_layout);
    control_Layout->addLayout(x_layout);
    control_Layout->addLayout(w_layout);
    control_Layout->addLayout(b_layout);
    control_Layout->addLayout(noise_layout);
    control_Layout->addWidget(m_startFit);

    // Create layout for grid and detached legend
    QHBoxLayout *mainLayout = new QHBoxLayout();
    mainLayout->addLayout(control_Layout);
    mainLayout->addWidget(m_chartView);
    mainLayout->setStretch(0, 0);
    mainLayout->setStretch(1, 10);
    setLayout(mainLayout);
}

void MainWidget::createSeries()
{
    // Fit series
    m_serieFit = new QLineSeries();
    m_serieFit->setColor(QColorConstants::Red);
    m_serieFit->setUseOpenGL(true);

    // data series
    m_serieData = new QScatterSeries();
    m_serieData->setMarkerShape(QScatterSeries::MarkerShapeCircle);
    m_serieData->setMarkerSize(4.0);
    m_serieData->setUseOpenGL(true);

    getSliderValue();

    // Create x/y vector for fit and start Gaussian
    N = 1000;
    datas = new double[N];

    for (int i = 0; i < N; i++) {
        // initialize data with noise
        datas[i] = param[0] * exp(-2 * pow(i - param[1], 2) / pow(param[2], 2)) + param[3];
        datas[i] += param[4] * param[0] * (rand() % 1000) / 1000.0;

        // add data to XY series
        m_serieData->append(i, datas[i]);
    }
    m_chart->addSeries(m_serieData);
    m_chart->addSeries(m_serieFit);
    work = new worker(N, 4, &mutex);
    connect(this, &MainWidget::seriesChange, work, &worker::setData, Qt::QueuedConnection);
    connect(work, &worker::fitEND, this, &MainWidget::copyParam, Qt::QueuedConnection);
}

void MainWidget::startFitting()
{
    work->startFitting(datas);
}

void MainWidget::getSliderValue()
{
    param[0] = (double) A_slider->value();
    param[1] = (double) x_slider->value();
    param[2] = (double) w_slider->value();
    param[3] = (double) b_slider->value();
    param[4] = noise_slider->value() / 100.0;

    /*
    std::clog << "Params " << std::endl;
    std::clog << "A      = " << param[0] << std::endl;
    std::clog << "x0     = " << param[1] << std::endl;
    std::clog << "w0     = " << param[2] << std::endl;
    std::clog << "b      = " << param[3] << std::endl;
*/
}

void MainWidget::copyParam(double *fit_params)
{
    std::scoped_lock lock(mutex);
    memcpy(fitParam, fit_params, sizeof(double) * 4);
    emit drawFit();
}

void MainWidget::drawingFit()
{
    m_serieFit->clear();
    double fit = 0.0;
    for (int i = 0; i <= 200; i++) {
        double I = i / 200.0 * N;
        // initialize data with noise
        fit = fitParam[0] * exp(-2 * pow(I - fitParam[1], 2) / pow(fitParam[2], 2)) + fitParam[3];

        // add data to XY series
        m_serieFit->append(I, fit);
    }
}

void MainWidget::updataSeries()
{
    std::scoped_lock lock(mutex);
    getSliderValue();
    m_serieData->clear();
    for (int i = 0; i < N; i++) {
        // initialize data with noise
        datas[i] = param[0] * exp(-2 * pow(i - param[1], 2) / pow(param[2], 2)) + param[3];
        datas[i] += param[4] * param[0] * (rand() % 1000) / 1000.0;

        // add data to XY series
        m_serieData->append(i, datas[i]);
    }
    emit seriesChange(datas);
}
