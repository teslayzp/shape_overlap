package knu.lsy.shapeoverlap.model;

import java.util.Objects;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonTypeName;

@JsonTypeName("regular_polygon")
public class RegularPolygon extends Shape {
    @JsonProperty("sides")
    private int sides;

    @JsonProperty("radius")
    private double radius;

    @JsonProperty("rotation")
    private double rotation;

    public RegularPolygon() {
        super();
    }

    public RegularPolygon(double x, double y, int sides, double radius, double rotation) {
        super(x, y);
        this.sides = sides;
        this.radius = radius;
        this.rotation = rotation;
    }

    public int getSides() {
        return sides;
    }

    public void setSides(int sides) {
        this.sides = sides;
    }

    public double getRadius() {
        return radius;
    }

    public void setRadius(double radius) {
        this.radius = radius;
    }

    public double getRotation() {
        return rotation;
    }

    public void setRotation(double rotation) {
        this.rotation = rotation;
    }

    /**
     * Get vertices of the regular polygon
     */
    public double[][] getVertices() {
        double[][] vertices = new double[sides][2];
        double angleStep = 2 * Math.PI / sides;

        for (int i = 0; i < sides; i++) {
            double angle = i * angleStep + Math.toRadians(rotation);
            vertices[i][0] = getX() + radius * Math.cos(angle);
            vertices[i][1] = getY() + radius * Math.sin(angle);
        }

        return vertices;
    }

    @Override
    public boolean overlaps(Shape other) {
        return switch (other) {
            case Circle circle -> circle.overlaps(this);
            case RegularPolygon regularPolygon -> overlapsWithPolygon(regularPolygon);
            case IrregularPolygon irregularPolygon -> overlapsWithIrregularPolygon(irregularPolygon);
            default -> false;
        };
    }

    /**
     * Check overlap with another regular polygon using SAT (Separating Axis Theorem)
     */
    private boolean overlapsWithPolygon(RegularPolygon other) {
        double[][] vertices1 = this.getVertices();
        double[][] vertices2 = other.getVertices();

        return satOverlapCheck(vertices1, vertices2);
    }

    /**
     * Check overlap with irregular polygon using SAT
     */
    private boolean overlapsWithIrregularPolygon(IrregularPolygon other) {
        double[][] vertices1 = this.getVertices();
        double[][] vertices2 = other.getVertices();

        return satOverlapCheck(vertices1, vertices2);
    }

    /**
     * SAT (Separating Axis Theorem) implementation
     * If we can find a separating axis, the polygons don't overlap
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
        double max = -Double.MAX_VALUE;  // Changed from Double.MIN_VALUE

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
        double[][] vertices = getVertices();
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
        // Area of regular polygon = (n * s^2) / (4 * tan(π/n))
        // where s is side length = 2 * radius * sin(π/n)
        double sideLength = 2 * radius * Math.sin(Math.PI / sides);
        return (sides * sideLength * sideLength) / (4 * Math.tan(Math.PI / sides));
    }

    @Override
    public String toString() {
        return String.format("RegularPolygon{x=%.2f, y=%.2f, sides=%d, radius=%.2f, rotation=%.2f}",
                getX(), getY(), sides, radius, rotation);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;
        RegularPolygon that = (RegularPolygon) o;
        return sides == that.sides &&
                Double.compare(that.radius, radius) == 0 &&
                Double.compare(that.rotation, rotation) == 0;
    }

    @Override
    public int hashCode() {
        return Objects.hash(super.hashCode(), sides, radius, rotation);
    }
}