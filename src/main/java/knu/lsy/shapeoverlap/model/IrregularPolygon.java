package knu.lsy.shapeoverlap.model;

import java.util.Arrays;
import java.util.Objects;

import com.fasterxml.jackson.annotation.JsonProperty;

public class IrregularPolygon extends Shape {
    @JsonProperty("vertices")
    private double[][] vertices;

    public IrregularPolygon() {
        super();
    }

    public IrregularPolygon(double x, double y, double[][] vertices) {
        super(x, y);
        this.vertices = translateVertices(vertices, x, y);
    }

    public double[][] getVertices() {
        return vertices;
    }

    public void setVertices(double[][] vertices) {
        this.vertices = vertices;
    }

    /**
     * Translate vertices relative to center position
     */
    private double[][] translateVertices(double[][] originalVertices, double centerX, double centerY) {
        double[][] translatedVertices = new double[originalVertices.length][2];
        for (int i = 0; i < originalVertices.length; i++) {
            translatedVertices[i][0] = originalVertices[i][0] + centerX;
            translatedVertices[i][1] = originalVertices[i][1] + centerY;
        }
        return translatedVertices;
    }

    @Override
    public boolean overlaps(Shape other) {
        return switch (other) {
            case Circle circle -> circle.overlaps(this);
            case RegularPolygon regularPolygon -> overlapsWithRegularPolygon(regularPolygon);
            case IrregularPolygon irregularPolygon -> overlapsWithIrregularPolygon(irregularPolygon);
            default -> false;
        };
    }

    /**
     * Check overlap with regular polygon using SAT
     */
    private boolean overlapsWithRegularPolygon(RegularPolygon other) {
        double[][] vertices1 = this.getVertices();
        double[][] vertices2 = other.getVertices();

        return satOverlapCheck(vertices1, vertices2);
    }

    /**
     * Check overlap with another irregular polygon using SAT
     */
    private boolean overlapsWithIrregularPolygon(IrregularPolygon other) {
        double[][] vertices1 = this.getVertices();
        double[][] vertices2 = other.getVertices();

        return satOverlapCheck(vertices1, vertices2);
    }

    /**
     * SAT (Separating Axis Theorem) implementation
     */
    private boolean satOverlapCheck(double[][] vertices1, double[][] vertices2) {
        // Check axes from first polygon
        if (!checkSeparatingAxes(vertices1, vertices1, vertices2)) {
            return false;
        }

        // Check axes from second polygon
        return checkSeparatingAxes(vertices2, vertices1, vertices2);
    }

    /**
     * Check if any axis from the given polygon separates the two polygons
     */
    private boolean checkSeparatingAxes(double[][] axesSource, double[][] poly1, double[][] poly2) {
        for (int i = 0; i < axesSource.length; i++) {
            int next = (i + 1) % axesSource.length;

            // Get perpendicular vector (normal) to the edge
            double edgeX = axesSource[next][0] - axesSource[i][0];
            double edgeY = axesSource[next][1] - axesSource[i][1];
            double normalX = -edgeY;
            double normalY = edgeX;

            // Normalize the normal vector
            double length = Math.sqrt(normalX * normalX + normalY * normalY);
            if (length > 0) {
                normalX /= length;
                normalY /= length;
            }

            // Project both polygons onto this axis
            double[] projection1 = projectPolygon(poly1, normalX, normalY);
            double[] projection2 = projectPolygon(poly2, normalX, normalY);

            // Check if projections overlap
            if (projection1[1] < projection2[0] || projection2[1] < projection1[0]) {
                return false; // Found separating axis
            }
        }
        return true; // No separating axis found on these axes
    }

    /**
     * Project polygon onto an axis and return [min, max] values
     */
    private double[] projectPolygon(double[][] vertices, double axisX, double axisY) {
        double min = Double.MAX_VALUE;
        double max = -Double.MAX_VALUE;

        for (double[] vertex : vertices) {
            double projection = vertex[0] * axisX + vertex[1] * axisY;
            min = Math.min(min, projection);
            max = Math.max(max, projection);
        }

        return new double[]{min, max};
    }

    /**
     * Check if a point is inside this polygon using ray casting algorithm
     */
    public boolean isPointInside(double px, double py) {
        if (vertices == null || vertices.length < 3) {
            return false;
        }

        boolean inside = false;

        for (int i = 0, j = vertices.length - 1; i < vertices.length; j = i++) {
            double xi = vertices[i][0], yi = vertices[i][1];
            double xj = vertices[j][0], yj = vertices[j][1];

            if (((yi > py) != (yj > py)) &&
                    (px < (xj - xi) * (py - yi) / (yj - yi) + xi)) {
                inside = !inside;
            }
        }

        return inside;
    }

    @Override
    public double getArea() {
        if (vertices == null || vertices.length < 3) {
            return 0.0;
        }

        // Shoelace formula for polygon area
        double area = 0.0;
        int n = vertices.length;

        for (int i = 0; i < n; i++) {
            int j = (i + 1) % n;
            area += vertices[i][0] * vertices[j][1];
            area -= vertices[j][0] * vertices[i][1];
        }

        return Math.abs(area) / 2.0;
    }

    /**
     * Get the centroid of the polygon
     */
    public double[] getCentroid() {
        if (vertices == null || vertices.length == 0) {
            return new double[]{getX(), getY()};
        }

        double cx = 0.0, cy = 0.0;
        for (double[] vertex : vertices) {
            cx += vertex[0];
            cy += vertex[1];
        }

        return new double[]{cx / vertices.length, cy / vertices.length};
    }

    /**
     * Check if this polygon is convex
     */
    public boolean isConvex() {
        if (vertices == null || vertices.length < 3) {
            return false;
        }

        boolean hasPositive = false;
        boolean hasNegative = false;
        int n = vertices.length;

        for (int i = 0; i < n; i++) {
            int prev = (i - 1 + n) % n;
            int next = (i + 1) % n;

            double dx1 = vertices[i][0] - vertices[prev][0];
            double dy1 = vertices[i][1] - vertices[prev][1];
            double dx2 = vertices[next][0] - vertices[i][0];
            double dy2 = vertices[next][1] - vertices[i][1];

            double crossProduct = dx1 * dy2 - dy1 * dx2;

            if (crossProduct > 0) {
                hasPositive = true;
            } else if (crossProduct < 0) {
                hasNegative = true;
            }

            if (hasPositive && hasNegative) {
                return false; // Found both positive and negative, so not convex
            }
        }

        return true;
    }

    @Override
    public String toString() {
        return String.format("IrregularPolygon{x=%.2f, y=%.2f, vertices=%s}",
                getX(), getY(), Arrays.deepToString(vertices));
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;
        IrregularPolygon that = (IrregularPolygon) o;
        return Arrays.deepEquals(vertices, that.vertices);
    }

    @Override
    public int hashCode() {
        int result = Objects.hash(super.hashCode());
        result = 31 * result + Arrays.deepHashCode(vertices);
        return result;
    }
}