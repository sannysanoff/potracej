import compat.ConvertToJavaCurves;
import compat.PathElement;
import potracej.*;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.plaf.metal.MetalButtonUI;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.GeneralPath;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: san
 * Date: 6/10/12
 * Time: 12:55 PM
 * To change this template use File | Settings | File Templates.
 */
public class Main {

    static BufferedImage result;
    static Bitmap bmp;
    static param_t param = new param_t();
    static double scale = 3;
    static ImageIcon resultIcon;
    static ImageIcon srcIcon;
    static BufferedImage sourceImage;
    static boolean renderSourceImage = false;

    public static void main(String[] args) throws IOException {

        sourceImage = ImageIO.read(new File("girl.png"));

        //Toolkit.getDefaultToolkit().
        WritableRaster raster = sourceImage.getRaster();
        int[] iarr = new int[4];
        bmp = new Bitmap((int)(sourceImage.getWidth()), (int)(sourceImage.getHeight()));
        for(int y=0; y<sourceImage.getHeight(); y++) {
            for(int x=0; x<sourceImage.getWidth(); x++) {
                int[] pixel = raster.getPixel(x, y, iarr);
                if (pixel[0] + pixel[1] + pixel[2] + pixel[3] != 0) {
                    bmp.put(x, y, 1);
                }
            }
        }

        BufferedImage d2 = new BufferedImage((int) (scale * sourceImage.getWidth()), (int)(scale * sourceImage.getHeight()), BufferedImage.TYPE_INT_ARGB);
        Graphics2D d2g = (Graphics2D) d2.getGraphics();
        d2g.scale(scale, scale);
        d2g.drawImage(sourceImage, 0, 0, null);
        d2g.dispose();
        sourceImage.flush();
        srcIcon = new ImageIcon(d2);

        doTrace(scale);


        JFrame frame = new JFrame("Result") {
            {
                setLayout(new BorderLayout());
                resultIcon = new ImageIcon(result, "Result");
                JButton resultButton = new JButton(resultIcon);
                resultButton.setUI(new MetalButtonUI() {
                    @Override
                    protected void paintButtonPressed(Graphics g, AbstractButton b) {
                        //
                    }
                });
                add(resultButton, BorderLayout.CENTER);
                resultButton.setPressedIcon(srcIcon);
                JPanel stuff = new JPanel();
                add(stuff, BorderLayout.NORTH);
                stuff.setLayout(new GridLayout(4, 2));
                stuff.add(new JLabel("Suppress speckles"));
                final JSlider turdSlider = new JSlider(JSlider.HORIZONTAL, 0, 100, param.turdsize);
                stuff.add(turdSlider);
                turdSlider.addChangeListener(new ChangeListener() {
                    @Override
                    public void stateChanged(ChangeEvent e) {
                        param.turdsize = turdSlider.getValue();
                        doRetrace();
                    }
                });
                stuff.add(new JLabel("Smooth corners"));
                final JSlider smoothSlider = new JSlider(JSlider.HORIZONTAL, 0, 300, (int) (param.opttolerance * 100));
                stuff.add(smoothSlider);
                smoothSlider.addChangeListener(new ChangeListener() {
                    @Override
                    public void stateChanged(ChangeEvent e) {
                        param.opttolerance = smoothSlider.getValue() / 100.0;
                        doRetrace();
                    }
                });
                stuff.add(new JLabel("Optimize paths"));
                final JSlider optSlider = new JSlider(JSlider.HORIZONTAL, 0, 125, (int) (param.alphamax * 100));
                stuff.add(optSlider);
                optSlider.addChangeListener(new ChangeListener() {
                    @Override
                    public void stateChanged(ChangeEvent e) {
                        param.alphamax = optSlider.getValue()/100.0;
                        doRetrace();
                    }
                });
                final JCheckBox renderSource = new JCheckBox("Render source");
                stuff.add(renderSource);
                renderSource.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        renderSourceImage = renderSource.getModel().isArmed();
                        doRetrace();
                    }
                });

            }

            private void doRetrace() {
                doTrace(scale);

                resultIcon.setImage(result);
                repaint();
            }
        };
        frame.pack();
        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
        frame.setVisible(true);

    }

    private static void doTrace(double scale) {
        PoTraceJ poTraceJ = new PoTraceJ(param);
        long l = System.currentTimeMillis();

        path_t trace = null;
        for(int i=0; i<10; i++) {
            trace = poTraceJ.trace(bmp);
            Thread.yield();
        }
        poTraceJ.resetTimers();
        for(int i=0; i<100; i++) {
            trace = poTraceJ.trace(bmp);
        }
        poTraceJ.printTimers();
        l = System.currentTimeMillis() - l;
        System.out.println("L="+l);
        ArrayList<PathElement> al = new ArrayList<PathElement>();
        ConvertToJavaCurves.convert(trace, new HashSet<ConvertToJavaCurves.Point>(), al);

        if (result != null)
            result.flush();
        result = new BufferedImage((int)(scale * bmp.getWidth()), (int)(scale * bmp.getHeight()), BufferedImage.TYPE_INT_ARGB);


        Graphics2D g2 = (Graphics2D)result.getGraphics();
        g2.scale(scale, scale);
        g2.setColor(Color.WHITE);
        g2.fillRect(0, 0, bmp.getWidth(), bmp.getHeight());
        g2.setColor(Color.BLACK);
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2.setRenderingHint(RenderingHints.KEY_FRACTIONALMETRICS, RenderingHints.VALUE_FRACTIONALMETRICS_ON);
        g2.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
        GeneralPath path = new GeneralPath();
        for (PathElement pathElement : al) {
            switch (pathElement.getType()) {
                case CLOSE_PATH:
                    path.closePath();
                    break;
                case LINE_TO:
                    path.lineTo(pathElement.getP0x(), pathElement.getP0y());
                    break;
                case MOVE_TO:
                    path.moveTo(pathElement.getP0x(), pathElement.getP0y());
                    break;
                case CURVE_TO:
                    path.curveTo(pathElement.getP0x(), pathElement.getP0y(), pathElement.getP1x(), pathElement.getP1y(), pathElement.getP2x(), pathElement.getP2y());
                    break;
            }
        }
        g2.setPaint(Color.black);
        g2.fill(path);
    }

}
