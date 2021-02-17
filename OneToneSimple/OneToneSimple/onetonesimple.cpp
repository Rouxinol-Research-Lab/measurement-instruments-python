#include "onetonesimple.h"

OneToneSimple::OneToneSimple(QWidget *parent)
    : QMainWindow(parent)
{
    ui.setupUi(this);

	ui.PlotWidget->yAxis->setRange(-110, 0);
	ui.PlotWidget->xAxis->setRange(4.9e9, 6.2e9);

	ui.PlotWidget->addGraph();
	ui.PlotWidget->graph(0)->setPen(QColor(0, 0, 255, 255));
	ui.PlotWidget->graph(0)->setLineStyle(QCPGraph::lsNone);
	ui.PlotWidget->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
	ui.PlotWidget->graph(0)->setName("Amplitude");
}

void OneToneSimple::receivedLog(char* log)
{
	qDebug("Thread id inside receivedLog %d", (int)QThread::currentThreadId());
	qDebug("message %s", log);
	ui.LogWidget->addItem(log);
}


void OneToneSimple::receivedDataPoint(double datax, double datay)
{
	qDebug("Thread id inside receivedDataPoint %d", (int)QThread::currentThreadId());
	qDebug("freq %f amp %f", datax, datay);
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
