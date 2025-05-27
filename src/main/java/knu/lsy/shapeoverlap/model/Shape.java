package knu.lsy.shapeoverlap.model;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonSubTypes;
import com.fasterxml.jackson.annotation.JsonTypeInfo;
import java.util.Objects;

@JsonTypeInfo(
        use = JsonTypeInfo.Id.NAME,
        include = JsonTypeInfo.As.PROPERTY,
        property = "type"
)
@JsonSubTypes({
        @JsonSubTypes.Type(value = Circle.class, name = "circle"),
        @JsonSubTypes.Type(value = RegularPolygon.class, name = "regular_polygon"),
        @JsonSubTypes.Type(value = IrregularPolygon.class, name = "irregular_polygon")
})
public abstract class Shape {
    @JsonProperty("x")
    private double x;

    @JsonProperty("y")
    private double y;

    public Shape() {
        this.x = 0.0;
        this.y = 0.0;
    }

    public Shape(double x, double y) {
        this.x = x;
        this.y = y;
    }

    public double getX() {
        return x;
    }

    public void setX(double x) {
        this.x = x;
    }

    public double getY() {
        return y;
    }

    public void setY(double y) {
        this.y = y;
    }

    /**
     * Abstract method to check if this shape overlaps with another shape
     * Each concrete shape class must implement this method
     */
    public abstract boolean overlaps(Shape other);

    /**
     * Abstract method to calculate the area of the shape
     */
    public abstract double getArea();

    /**
     * Calculate distance between this shape's center and another shape's center
     */
    public double distanceTo(Shape other) {
        double dx = this.x - other.x;
        double dy = this.y - other.y;
        return Math.sqrt(dx * dx + dy * dy);
    }

    /**
     * Move the shape by given offsets
     */
    public void translate(double dx, double dy) {
        this.x += dx;
        this.y += dy;
    }

    /**
     * Get the type of shape as string
     */
    public String getShapeType() {
        return this.getClass().getSimpleName();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Shape shape = (Shape) o;
        return Double.compare(shape.x, x) == 0 &&
                Double.compare(shape.y, y) == 0;
    }

    @Override
    public int hashCode() {
        return Objects.hash(x, y);
    }

    @Override
    public String toString() {
        return String.format("Shape{x=%.2f, y=%.2f}", x, y);
    }
}