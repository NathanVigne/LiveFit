#ifndef WORKER_H
#define WORKER_H

#include <QThread>
#include <QTimer>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>
#include <iostream>
#include <mutex>

struct data
{
    size_t n;
    double *t;
    double *y;
};

class worker : public QObject
{
    Q_OBJECT;

public:
    worker(const int N, const int P, std::mutex *mutex, QObject *parent = nullptr);
    ~worker();

    void startFitting(double *datas);
    void setMutex(std::mutex *newMutex);
    int setData(double *datas);

signals:
    void fitEND(double *params);

private slots:
    void copyData(double *datas);

private:
    void Fitting();
    void saveFitParam();
    int initialGuess();

    static int expb_f(const gsl_vector *x, void *data, gsl_vector *f);
    static int expb_df(const gsl_vector *x, void *data, gsl_matrix *f);

    const gsl_multifit_nlinear_type *T;
    gsl_multifit_nlinear_parameters params;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_fdf fdf;

    const size_t n;
    const size_t p;
    double *m_ydata;
    double *m_xdata;
    double *m_fitParam;
    data d;
    double x0, w0, A, b;
    double *param_guess;
    gsl_vector_view init_p;
    gsl_vector *f;

    double chisq, chisq0;
    int status, info;
    size_t i;

    QThread thread;
    std::mutex *m_mutex;
    QTimer *timer;
};

#endif // WORKER_H
