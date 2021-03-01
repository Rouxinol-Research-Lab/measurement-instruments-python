#include "onetonesimple.h"

OneToneSimple::OneToneSimple(QWidget *parent)
    : QMainWindow(parent)
{
    ui.setupUi(this);

	ui.PlotWidget->yAxis->setRange(-200, 0);

	ui.PlotWidget->addGraph();
	ui.PlotWidget->graph(0)->setPen(QColor(0, 0, 255, 255));
	ui.PlotWidget->graph(0)->setLineStyle(QCPGraph::lsNone);
	ui.PlotWidget->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
	ui.PlotWidget->graph(0)->setName("Amplitude");
}

void OneToneSimple::defineRange(double l, double r)
{
	ui.PlotWidget->xAxis->setRange(l, r);
	ui.PlotWidget->replot();
}

void OneToneSimple::receivedLog(const char* log)
{
	qDebug("Thread id inside receivedLog %d", (int)QThread::currentThreadId());
	ui.LogWidget->addItem(log);
}


void OneToneSimple::receivedDataPoint(double datax, double datay)
{
	qDebug("Thread id inside receivedDataPoint %d", (int)QThread::currentThreadId());
	qDebug("x: %f y %f", datax,datay);
	x.push_back(datax);
	y.push_back(datay);

	ui.PlotWidget->graph(0)->setData(x, y);

	ui.PlotWidget->replot();
}

void OneToneSimple::on_StopButton_clicked()
{
	qDebug("Thread id inside on_StopMeasurementButton_clicked %d", (int)QThread::currentThreadId());
	emit stopMeasurement();

}
