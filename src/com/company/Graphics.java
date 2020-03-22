package com.company;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYSplineRenderer;
import org.jfree.chart.title.TextTitle;
import org.jfree.chart.ui.ApplicationFrame;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.chart.ui.RectangleInsets;
import org.jfree.chart.ui.VerticalAlignment;
import org.jfree.chart.util.UnitType;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import java.awt.*;
import java.io.Console;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

//Базовий клас для створення main frame для простих додатків
public class Graphics extends ApplicationFrame {
    //Універсальний ідентифікатор версії для Serializable класу.
    //Десеріалізація використовує це число для того, щоб завантажений клас точно відповідав серіалізованому об'єкту.
    private static final long serialVersionUID = 3L;

    int N = 512;
    double[] result = new double[N];
    double[] T1 = new double[N];
    double[] T2 = new double[N];

    public Graphics(final String title) {
        super(title);

        JFreeChart chart = createChart();
        //Класс ChartPanel используется в качестве компонента графического интерфейса для отображения объекта JfreeChart
        //Этот конструктор создает панель, которая отображает указанную диаграмму.
        ChartPanel chartPanel = new ChartPanel(chart);

        chartPanel.setPreferredSize(new java.awt.Dimension(900, 700));
        // це дозволяє JToolBar навести курсор миші на ContentPanel
        setContentPane(chartPanel);
    }

    public Graphics(final String title, double[] res, int idx) {
        super(title);

        JFreeChart chart;
//        if (idx == 0) chart = createChart2(res);
//        else if (idx == 2) chart = createChart3();
//        else  chart = createChart4(res);
        chart = createChart5();
        ChartPanel chartPanel = new ChartPanel(chart);

        chartPanel.setPreferredSize(new java.awt.Dimension(900, 700));
        setContentPane(chartPanel);
    }

    private JFreeChart createChart() {
        /*Створює лінійну діаграму
        title- назва діаграми ( nullдозволено).
        xAxisLabel- мітка для осі X ( nullдозволено).
        yAxisLabel- мітка для осі Y ( nullдозволено).
        dataset- набір даних для діаграми ( nullдозволено).
        orientation- орієнтація ділянки (горизонтальна або вертикальна) ( nullНЕ дозволено).
        legend - прапор із зазначенням того, потрібна чи ні легенда.
        tooltips - налаштувати діаграму для генерації підказок інструментів?
        urls - налаштувати діаграму для генерування URL-адрес?*/
        final JFreeChart chart = ChartFactory.createXYLineChart(
                "Графік випадкового сигналу",
                null,
                null,
                null,
                PlotOrientation.VERTICAL,
                true,
                false,
                false
        );

        chart.setBackgroundPaint(new Color(153, 255, 102));
        //Загальний клас для побудови графіків даних у вигляді (x, y) пар.
        //Повертає plot, поданий у вигляді XYPlot.
        final XYPlot plot = chart.getXYPlot();
        
        plot.setBackgroundPaint(new Color(204, 255, 204));
        //Встановлює фарбу для ліній сітки, намічених проти осі домену
        plot.setDomainGridlinePaint(new Color(0, 51, 0));
        plot.setRangeGridlinePaint(new Color(0, 51, 0));
        //Встановлює зміщення осей прямокутника(графіка)
        plot.setAxisOffset(new RectangleInsets(1.0, 1.0, 1.0, 1.0));
        //Повертає вісь домену з індексом 0.
        ValueAxis axis = plot.getDomainAxis();
        //Встановлює прапор, який контролює, чи видно лінію осі чи ні
        axis.setAxisLineVisible(false);
        //Этот класс может иметь доступ к числовым данным любой оси.
        //getRangeAxis() Повертає вісь діапазону для діаграми.
        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        rangeAxis.setAxisLineVisible(false);
        //становлює джерело для отримання стандартних TickUnits для осі
        //Вісь намагатиметься вибрати найменший відміток від джерела, який не спричиняє перекриття міток
        rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        //Візуалізатор, який з'єднує точки даних із природними кубічними сплайнами та малює фігури в кожній точці даних
        XYSplineRenderer r0 = new XYSplineRenderer();
        r0.setPrecision(1);
        //Встановлює прапор "форми видимі" для точок даних
        r0.setSeriesShapesVisible(0, false);

        int n = 14, w = 2000;
        double A, q, temp_sum, Mx, Dx;

        for (int t = 0; t < N; t++) {
            temp_sum = 0;
            A = Math.random();
            q = Math.random();
            for (int k = 1; k <= n; k++) {
                temp_sum += A * Math.sin((w / n) * k * t + q);
            }
            result[t] = temp_sum;
        }

        XYSeries series1 = new XYSeries("x(t)");
        for (int i = 0; i < N; i++) {
            series1.add((double) i, result[i]);
        }

        temp_sum = 0;
        for (int i = 0; i < N; i++) {
            temp_sum += result[i];
        }
        Mx = temp_sum / N;

        temp_sum = 0;
        for (int i = 0; i < N; i++) {
            temp_sum += (Mx - result[i]) * (Mx - result[i]);
        }
        Dx = temp_sum / (N - 1);

        //Это класс, который представляет последовательность из нуля или более элементов данных в форме (x, y).
        //Этот конструктор создает новую пустую серию.
        XYSeries series2 = new XYSeries("Mx");
        series2.add(0, Mx);
        series2.add(N, Mx);

        XYSeries series3 = new XYSeries("Dx");
        series3.add(0, Dx);
        series3.add(N, Dx);

        TextTitle localTextTitle = new TextTitle("Mx = " + Mx);
        localTextTitle.setPosition(RectangleEdge.BOTTOM);
        localTextTitle.setPadding(new RectangleInsets(UnitType.RELATIVE, 0.05D, 0.05D, 0.05D, 0.05D));
        localTextTitle.setVerticalAlignment(VerticalAlignment.BOTTOM);
        chart.addSubtitle(localTextTitle);

        TextTitle localTextTitle2 = new TextTitle("Dx = " + Dx);
        localTextTitle2.setPosition(RectangleEdge.BOTTOM);
        localTextTitle2.setPadding(new RectangleInsets(UnitType.RELATIVE, 0.05D, 0.05D, 0.05D, 0.05D));
        localTextTitle2.setVerticalAlignment(VerticalAlignment.BOTTOM);
        chart.addSubtitle(localTextTitle2);

        //Представляє колекцію XYSeries об'єктів, які можна використовувати як набір даних.
        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series1);
        dataset.addSeries(series2);
        dataset.addSeries(series3);

        //Встановлює набір даних для сюжету та надсилає подію зміни всім зареєстрованим слухачам.
        plot.setDataset(0, dataset);

        //Встановлює візуалізатор для набору даних із вказаним індексом та надсилає подію зміни всім зареєстрованим слухачам
        plot.setRenderer(0, r0);

        //Встановлює фарбу, використовувану для серії
        r0.setSeriesPaint(0, Color.black);
        r0.setSeriesPaint(1, Color.red);
        r0.setSeriesPaint(2, Color.orange);

        return chart;
    }

    private double[] createChart2(double[] res) {
//        final JFreeChart chart = ChartFactory.createXYLineChart(
//                "Функція ДПФ",
//                null,
//                null,
//                null,
//                PlotOrientation.VERTICAL,
//                true,
//                false,
//                false
//        );
//
//        chart.setBackgroundPaint(new Color(153, 255, 102));
//
//        final XYPlot plot = chart.getXYPlot();
//        plot.setBackgroundPaint(new Color(204, 255, 204));
//
//        plot.setDomainGridlinePaint(new Color(0, 51, 0));
//        plot.setRangeGridlinePaint(new Color(0, 51, 0));
//
//        plot.setAxisOffset(new RectangleInsets(1.0, 1.0, 1.0, 1.0));
//
//        ValueAxis axis = plot.getDomainAxis();
//        axis.setAxisLineVisible(false);
//
//        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
//        rangeAxis.setAxisLineVisible(false);
//        rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
//
//        XYSplineRenderer r0 = new XYSplineRenderer();
//        r0.setPrecision(8);
//        r0.setSeriesShapesVisible(0, false);

        double F;
        double F_Re;
        double F_Im;
        double F_Im1, F_Im2;
        double F_Re1, F_Re2;

        XYSeries series1 = new XYSeries("T1");
        XYSeries series2 = new XYSeries("T2");

        double start;
        double end;
        double start2;
        double end2;
        double start3;
        double end3;
        double[] Mcl = new double[N];
        double[] Mt = new double[N];
        double[] Mf = new double[N];

        double[] _time = new double[N];
        double _start;
        double _end;
        for (int NN = 0; NN < N; NN++) {
            _start = System.nanoTime();
            start = System.nanoTime();
            for (int p = 0; p < NN; p++) {
                F_Re = 0;
                F_Im = 0;
                for (int k = 0; k < (NN - 1); k++) {
                    F_Re += res[k] * Math.cos(2 * Math.PI * p * k / NN);
                    F_Im += res[k] * Math.sin(2 * Math.PI * p * k / NN);
                }
                F = Math.sqrt((F_Re * F_Re) + (F_Im * F_Im));
            }
            end = System.nanoTime();
            Mcl[NN] = end - start;

            double[][] Wreal = new double[NN][NN];
            double[][] Wim = new double[NN][NN];

            for (int i = 0; i < NN; i++) {
                for (int j = 0; j < (NN - 1); j++) {
                    Wreal[i][j] = Math.cos(2 * Math.PI * i * j / NN);
                    Wim[i][j] = Math.sin(2 * Math.PI * i * j / NN);
                }
            }

            start2 = System.nanoTime();
            for (int p = 0; p < NN; p++) {
                F_Re = 0;
                F_Im = 0;
                for (int k = 0; k < (NN - 1); k++) {
                    F_Re += res[k] * Wreal[p][k];
                    F_Im += res[k] * Wim[p][k];
                }
                F = Math.sqrt((F_Re * F_Re) + (F_Im * F_Im));
            }
            end2 = System.nanoTime();
            Mt[NN] = end2 - start2;

            start3 = System.nanoTime();
            for (int p = 0; p < NN; p++) {
                F_Re1 = 0;
                F_Im1 = 0;
                F_Re2 = 0;
                F_Im2 = 0;
                for (int k = 0; k < (NN / 2) - 1; k++) {
                    F_Re2 += res[2 * k] * Math.cos(2 * Math.PI * p * k / (NN / 2));
                    F_Im2 += res[2 * k] * Math.sin(2 * Math.PI * p * k / (NN / 2));

                    F_Re1 += res[2 * k + 1] * Math.cos(2 * Math.PI * p * (2 * k + 1) / NN);
                    F_Im1 += res[2 * k + 1] * Math.sin(2 * Math.PI * p * (2 * k + 1) / NN);
                }
                F_Re = F_Re2 + F_Re1;
                F_Im = F_Im2 + F_Im1;
                F = Math.sqrt((F_Re * F_Re) + (F_Im * F_Im));
            }
            end3 = System.nanoTime();
            Mf[NN] = end3 - start3;
            _end = System.nanoTime();
            _time[NN] = _end - _start;
        }

//        for (int i = 0; i < N; i++) {
//            series1.add(i, Mcl[i] / Mt[i]);
//            series2.add(i, Mcl[i] / Mf[i]);
//            System.out.println(_time[i]);
//        }
//
//
//
//        //for (int p = 0; p < N; p++) {
//        //    series2.add(p, F[p]);
//        //}
//
//        XYSeriesCollection dataset = new XYSeriesCollection();
//        dataset.addSeries(series1);
//        dataset.addSeries(series2);
//
//        plot.setDataset(0, dataset);
//
//        plot.setRenderer(0, r0);
//
//        r0.setSeriesPaint(0, Color.black);

        return _time;
    }

    private JFreeChart createChart3() {
        final JFreeChart chart = ChartFactory.createXYLineChart(
                "Часова складність",
                null,
                null,
                null,
                PlotOrientation.VERTICAL,
                true,
                false,
                false
        );

        chart.setBackgroundPaint(new Color(153, 255, 102));

        final XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(new Color(204, 255, 204));

        plot.setDomainGridlinePaint(new Color(0, 51, 0));
        plot.setRangeGridlinePaint(new Color(0, 51, 0));

        plot.setAxisOffset(new RectangleInsets(1.0, 1.0, 1.0, 1.0));

        ValueAxis axis = plot.getDomainAxis();
        axis.setAxisLineVisible(false);

        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        rangeAxis.setAxisLineVisible(false);
        rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());

        XYSplineRenderer r0 = new XYSplineRenderer();
        r0.setPrecision(8);
        r0.setSeriesShapesVisible(0, false);

        ArrayList<Double> F = new ArrayList<Double>(10 * N);
        ArrayList<Double> F_Re = new ArrayList<Double>(10 * N);
        ArrayList<Double> F_Im = new ArrayList<Double>(10 * N);
        for (int i = 0; i < 10 * N; i++) {
            F.add(0.0); F_Re.add(0.0); F_Im.add(0.0);
        }

//        double[] F = new double[10*N];
//        double[] F_Re = new double[10*N];
//        double[] F_Im= new double[10*N];
        int n = 14, w = 2000;
        double A, q, temp_sum, start, end;

        XYSeries series3 = new XYSeries("t(N)");

        for (int N_new = 0; N_new < 2000;) {
            double[] res = new double[N_new];

            for (int t = 0; t < N_new; t++) {
                temp_sum = 0;
                A = Math.random();
                q = Math.random();
                for (int k = 1; k <= n; k++) {
                    temp_sum += A * Math.sin((w / n) * k * t + q);
                }
                res[t] = temp_sum;
            }

            start = System.nanoTime();
            for (int p = 0; p < N_new; p++) {
                for (int k = 0; k < (N_new - 1); k++) {
                    double fi = 2 * Math.PI * p * k / N_new;
                    F_Re.set(p, F_Re.get(p) + res[k] * Math.cos(fi));
                    F_Im.set((p), F_Im.get(p) + res[k] * Math.sin(fi));
                }
                F.set(p, Math.sqrt(Math.pow(F_Re.get(p), 2) + Math.pow(F_Im.get(p), 2)));
            }
            end = System.nanoTime();
            series3.add(N_new, end - start);
            //System.out.println(N_new);
            N_new += 20;
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        //dataset.addSeries(series1);
        dataset.addSeries(series3);

        plot.setDataset(0, dataset);

        plot.setRenderer(0, r0);

        r0.setSeriesPaint(0, Color.black);
        r0.setSeriesPaint(1, Color.blue);

        return chart;
    }

    private  JFreeChart createChart5() {
        final JFreeChart chart = ChartFactory.createXYLineChart(
                "Функція ДПФ / Функція  ШПФ",
                null,
                null,
                null,
                PlotOrientation.VERTICAL,
                true,
                false,
                false
        );

        chart.setBackgroundPaint(new Color(153, 255, 102));

        final XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(new Color(204, 255, 204));

        plot.setDomainGridlinePaint(new Color(0, 51, 0));
        plot.setRangeGridlinePaint(new Color(0, 51, 0));

        plot.setAxisOffset(new RectangleInsets(1.0, 1.0, 1.0, 1.0));

        ValueAxis axis = plot.getDomainAxis();
        axis.setAxisLineVisible(false);

        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        rangeAxis.setAxisLineVisible(false);
        //rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());

        XYSplineRenderer r0 = new XYSplineRenderer();
        r0.setPrecision(8);
        r0.setSeriesShapesVisible(0, false);

        XYSeries series4 = new XYSeries("T(N)");
        double[] T1 = createChart2(new double[N]);
        double[] T2 = createChart4(new double[N]);
        double[] T = new double[N];
        for (int p = 0; p < N; p++) {
            T[p] = T1[p] / T2[p];
            series4.add(p, T[p]);
        }
        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(series4);

        plot.setDataset(0, dataset);

        plot.setRenderer(0, r0);

        r0.setSeriesPaint(0, Color.black);
        r0.setSeriesPaint(1, Color.blue);

        return chart;
    }

    private double[] createChart4(double[] res) {
//        final JFreeChart chart = ChartFactory.createXYLineChart(
//                "Функція ШПФ",
//                null,
//                null,
//                null,
//                PlotOrientation.VERTICAL,
//                true,
//                false,
//                false
//        );
//
//        chart.setBackgroundPaint(new Color(153, 255, 102));
//
//        final XYPlot plot = chart.getXYPlot();
//        plot.setBackgroundPaint(new Color(204, 255, 204));
//
//        plot.setDomainGridlinePaint(new Color(0, 51, 0));
//        plot.setRangeGridlinePaint(new Color(0, 51, 0));
//
//        plot.setAxisOffset(new RectangleInsets(1.0, 1.0, 1.0, 1.0));
//
//        ValueAxis axis = plot.getDomainAxis();
//        axis.setAxisLineVisible(false);
//
//        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
//        rangeAxis.setAxisLineVisible(false);
//        //rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
//
//        XYSplineRenderer r0 = new XYSplineRenderer();
//        r0.setPrecision(8);
//        r0.setSeriesShapesVisible(0, false);
//
//        //double[] F = new double[N];

        XYSeries series4 = new XYSeries("F(p)");

        double[] time = new double[N];
        double start;
        double end;
        double[] F_Re2 = new double[N];
        double[] F_Im2 = new double[N];
        double[] F_Re1 = new double[N];
        double[] F_Im1 = new double[N];
        double[] F_Re = new double[N];
        double[] F_Im = new double[N];
        double[] F = new double[N];

        for (int p = 0; p < N; p++) {
            start = System.nanoTime();
            for (int k = 0; k < N/2 - 1; k++) {
                double fi = (2 * Math.PI * p * k) / (N/2);
                double fi2 = (2 * Math.PI * p * (2*k+1)) / (N);
                F_Re2[p] += res[2*k] * Math.cos(fi);
                F_Im2[p] += res[2*k] * Math.sin(fi);
                F_Re1[p] += res[2*k+1] * Math.cos(fi2);
                F_Im1[p] += res[2*k+1] * Math.sin(fi2);
            }
            F_Re[p] = F_Re2[p] + F_Re1[p];
            F_Im[p] = F_Im2[p] + F_Im1[p];
            F[p] = Math.sqrt(Math.pow(F_Re[p], 2) + Math.pow(F_Im[p], 2));
            end = System.nanoTime();
            time[p] = end - start;
        }
//        for (int p = 0; p < N; p++) {
//            series4.add(p, F[p]);
//            System.out.println(time[p]);
//        }
//        XYSeriesCollection dataset = new XYSeriesCollection();
//        //dataset.addSeries(series1);
//        dataset.addSeries(series4);
//
//        plot.setDataset(0, dataset);
//
//        plot.setRenderer(0, r0);
//
//        r0.setSeriesPaint(0, Color.black);
//        r0.setSeriesPaint(1, Color.blue);

        return time;
    }

}