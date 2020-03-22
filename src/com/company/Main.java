package com.company;

import org.jfree.ui.RefineryUtilities;

public class Main {

    public static void main(String[] args) {
        final Graphics demo = new Graphics("Функція ДПФ до шпф");
        demo.pack();
        RefineryUtilities.positionFrameRandomly(demo);
        demo.setVisible(true);
        double[] res = demo.result;
//
        final Graphics demo2 = new Graphics("Функція ДПФ до шпф", res, 0);
        demo2.pack();
        RefineryUtilities.positionFrameRandomly(demo2);
        demo2.setVisible(true);
//
//        final Graphics demo3 = new Graphics("Графік ШПФ", res, 1);
//        demo3.pack();
//        RefineryUtilities.positionFrameRandomly(demo3);
//        demo3.setVisible(true);
//
//        final Graphics demo4 = new Graphics("Графік часової складності", res, 2);
//        demo4.pack();
//        RefineryUtilities.positionFrameRandomly(demo4);
//        demo4.setVisible(true);

    }
}