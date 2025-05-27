package knu.lsy.shapeoverlap.model;

import java.util.Objects;

import com.fasterxml.jackson.annotation.JsonProperty;

public class Circle extends Shape {
    @JsonProperty("radius")
    private double radius;

    public Circle() {
        super();
    }

    public Circle(double x, double y, double radius) {
        super(x, y);
        this.radius = radius;
    }

    public double getRadius() {
        return radius;
    }

    public void setRadius(double radius) {
        this.radius = radius;
    }

    @Override
    public boolean overlaps(Shape other) {
        if (other instanceof Circle circle) {
            return overlapsWithCircle(circle);
        } else if (other instanceof RegularPolygon || other instanceof IrregularPolygon) {
            return overlapsWithPolygon(other);
        }
        return false;
    }

    /**
     * Check if this circle overlaps with another circle
     * Formula: distance between centers < sum of radii
     */
    private boolean overlapsWithCircle(Circle other) {
        double dx = this.getX() - other.getX();
        double dy = this.getY() - other.getY();
        double distance = Math.sqrt(dx * dx + dy * dy);
        return distance < (this.radius + other.getRadius());
    }

    /**
     * Check if this circle overlaps with any polygon type
     * Uses combination of vertex-in-circle and edge-circle intersection tests
     */
    private boolean overlapsWithPolygon(Shape polygon) {
        if (polygon == null) {
            throw new IllegalArgumentException("Polygon cannot be null");
        }
        
        double[][] vertices;
        switch (polygon) {
            case RegularPolygon regularPolygon -> vertices = regularPolygon.getVertices();
            case IrregularPolygon irregularPolygon -> vertices = irregularPolygon.getVertices();
            default -> throw new IllegalArgumentException("Unsupported polygon type: " + polygon.getClass().getName());
        }

        // Check if any vertex of the polygon is inside the circle
        for (double[] vertex : vertices) {
            if (isPointInCircle(vertex[0], vertex[1])) {
                return true;
            }
        }

        // Check if circle center is inside the polygon
        boolean isInside = switch (polygon) {
            case RegularPolygon regularPolygon -> regularPolygon.isPointInside(this.getX(), this.getY());
            case IrregularPolygon irregularPolygon -> irregularPolygon.isPointInside(this.getX(), this.getY());
            default -> throw new IllegalArgumentException("Unsupported polygon type: " + polygon.getClass().getName());
        };
        if (isInside) {
            return true;
        }

        // Check if any edge of the polygon intersects with the circle
        for (int i = 0; i < vertices.length; i++) {
            int next = (i + 1) % vertices.length;
            if (circleIntersectsLineSegment(
                    vertices[i][0], vertices[i][1],
                    vertices[next][0], vertices[next][1]
            )) {
                return true;
            }
        }

        return false;
    }

    /**
     * Check if a point is inside this circle
     */
    private boolean isPointInCircle(double px, double py) {
        double dx = px - this.getX();
        double dy = py - this.getY();
        return (dx * dx + dy * dy) <= (radius * radius);
    }

    /**
     * Check if circle intersects with a line segment
     * Uses distance from circle center to line segment
     */
    private boolean circleIntersectsLineSegment(double x1, double y1, double x2, double y2) {
        double cx = this.getX();
        double cy = this.getY();

        // Vector from point 1 to point 2
        double dx = x2 - x1;
        double dy = y2 - y1;

        // Vector from point 1 to circle center
        double fx = cx - x1;
        double fy = cy - y1;

        double a = dx * dx + dy * dy;
        double b = 2 * (fx * dx + fy * dy);
        double c = (fx * fx + fy * fy) - radius * radius;

        double discriminant = b * b - 4 * a * c;

        if (discriminant < 0) {
            return false; // No intersection
        }

        discriminant = Math.sqrt(discriminant);

        double t1 = (-b - discriminant) / (2 * a);
        double t2 = (-b + discriminant) / (2 * a);

        // Check if intersection points are within the line segment
        return (t1 >= 0 && t1 <= 1) || (t2 >= 0 && t2 <= 1) || (t1 < 0 && t2 > 1);
    }

    @Override
    public double getArea() {
        return Math.PI * radius * radius;
    }

    @Override
    public String toString() {
        return String.format("Circle{x=%.2f, y=%.2f, radius=%.2f}", getX(), getY(), radius);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        if (!super.equals(o)) return false;
        Circle circle = (Circle) o;
        return Double.compare(circle.radius, radius) == 0;
    }

    @Override
    public int hashCode() {
        return Objects.hash(super.hashCode(), radius);
    }
}