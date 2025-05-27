
package knu.lsy.shapeoverlap.validation;

import knu.lsy.shapeoverlap.model.Circle;
import knu.lsy.shapeoverlap.model.RegularPolygon;
import knu.lsy.shapeoverlap.model.IrregularPolygon;
import knu.lsy.shapeoverlap.exception.ValidationException;

public class ShapeValidator {

    public static void validateCircle(Circle circle) {
        if (circle.getRadius() <= 0) {
            throw new ValidationException("Circle radius must be positive");
        }
    }

    public static void validateRegularPolygon(RegularPolygon polygon) {
        if (polygon.getSides() < 3) {
            throw new ValidationException("Regular polygon must have at least 3 sides");
        }
        if (polygon.getRadius() <= 0) {
            throw new ValidationException("Regular polygon radius must be positive");
        }
    }

    public static void validateIrregularPolygon(IrregularPolygon polygon) {
        if (polygon.getVertices() == null || polygon.getVertices().length < 3) {
            throw new ValidationException("Irregular polygon must have at least 3 vertices");
        }
        // Check for valid vertex coordinates
        for (double[] vertex : polygon.getVertices()) {
            if (vertex == null || vertex.length != 2) {
                throw new ValidationException("Invalid vertex coordinates");
            }
        }
    }
}