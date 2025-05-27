package knu.lsy.shapeoverlap.controller;

import knu.lsy.shapeoverlap.model.Shape;
import knu.lsy.shapeoverlap.service.OverlapService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.*;
import io.swagger.v3.oas.annotations.Operation;
import io.swagger.v3.oas.annotations.tags.Tag;
import jakarta.validation.Valid;
import jakarta.validation.constraints.NotEmpty;
import jakarta.validation.constraints.Min;

import java.util.List;
import java.util.Map;
import java.util.HashMap;

@RestController
@RequestMapping("/api/shapes")
@CrossOrigin(origins = "*")
@Tag(name = "Shape Operations", description = "API endpoints for shape overlap detection and analysis")
public class ShapeController {

    @Autowired
    private OverlapService overlapService;

    @Operation(summary = "Check overlaps between shapes", 
              description = "Analyzes a list of shapes and returns detailed overlap information")
    @PostMapping("/check-overlaps")
    public ResponseEntity<Map<String, Object>> checkOverlaps(@Valid @NotEmpty @RequestBody List<Shape> shapes) {
        try {
            // Find overlaps
            List<OverlapService.OverlapResult> overlaps = overlapService.findOverlaps(shapes);

            // Group overlapping shapes
            List<List<Integer>> groups = overlapService.groupOverlappingShapes(shapes);

            // Calculate statistics
            OverlapService.OverlapStatistics statistics = overlapService.calculateStatistics(shapes);

            // Prepare response
            Map<String, Object> response = new HashMap<>();
            response.put("success", true);
            response.put("overlaps", overlaps);
            response.put("groups", groups);
            response.put("statistics", statistics);
            response.put("totalShapes", shapes.size());
            response.put("overlappingPairs", overlaps.size());

            return ResponseEntity.ok(response);

        } catch (Exception e) {
            return ResponseEntity.internalServerError()
                    .body(createErrorResponse("Error processing shapes: " + e.getMessage()));
        }
    }

    @Operation(summary = "Count overlaps between shapes",
              description = "Returns only the count of overlapping shape pairs")
    @PostMapping("/count-overlaps")
    public ResponseEntity<Map<String, Object>> countOverlaps(@Valid @NotEmpty @RequestBody List<Shape> shapes) {
        try {
            int overlapCount = overlapService.countOverlaps(shapes);

            Map<String, Object> response = new HashMap<>();
            response.put("success", true);
            response.put("totalShapes", shapes.size());
            response.put("overlappingPairs", overlapCount);
            response.put("hasOverlaps", overlapCount > 0);

            return ResponseEntity.ok(response);

        } catch (Exception e) {
            return ResponseEntity.internalServerError()
                    .body(createErrorResponse("Error counting overlaps: " + e.getMessage()));
        }
    }

    @Operation(summary = "Find shapes overlapping with target",
              description = "Finds all shapes that overlap with a specific target shape")
    @PostMapping("/find-overlapping/{targetIndex}")
    public ResponseEntity<Map<String, Object>> findOverlappingShapes(
            @PathVariable @Min(0) int targetIndex,
            @Valid @NotEmpty @RequestBody List<Shape> shapes) {

        try {
            if (targetIndex >= shapes.size()) {
                return ResponseEntity.badRequest()
                        .body(createErrorResponse("Invalid target shape index"));
            }

            Shape targetShape = shapes.get(targetIndex);
            List<Shape> overlappingShapes = overlapService.findOverlappingShapes(targetShape, shapes);

            Map<String, Object> response = new HashMap<>();
            response.put("success", true);
            response.put("targetShape", targetShape);
            response.put("targetIndex", targetIndex);
            response.put("overlappingShapes", overlappingShapes);
            response.put("overlappingCount", overlappingShapes.size());

            return ResponseEntity.ok(response);

        } catch (Exception e) {
            return ResponseEntity.internalServerError()
                    .body(createErrorResponse("Error finding overlapping shapes: " + e.getMessage()));
        }
    }

    @Operation(summary = "Get overlap statistics",
              description = "Calculates detailed statistics about shape overlaps")
    @PostMapping("/statistics")
    public ResponseEntity<Map<String, Object>> getStatistics(@Valid @NotEmpty @RequestBody List<Shape> shapes) {
        try {
            OverlapService.OverlapStatistics statistics = overlapService.calculateStatistics(shapes);

            Map<String, Object> response = new HashMap<>();
            response.put("success", true);
            response.put("statistics", statistics);

            return ResponseEntity.ok(response);

        } catch (Exception e) {
            return ResponseEntity.internalServerError()
                    .body(createErrorResponse("Error calculating statistics: " + e.getMessage()));
        }
    }

    @Operation(summary = "Validate shapes",
              description = "Validates the input shapes and returns information about each shape")
    @PostMapping("/validate")
    public ResponseEntity<Map<String, Object>> validateShapes(@Valid @RequestBody List<Shape> shapes) {
        try {
            List<Map<String, Object>> shapeInfo = shapes.stream()
                    .map(this::getShapeInfo)
                    .toList();

            Map<String, Object> response = new HashMap<>();
            response.put("success", true);
            response.put("totalShapes", shapes.size());
            response.put("valid", true);
            response.put("shapes", shapeInfo);

            return ResponseEntity.ok(response);

        } catch (Exception e) {
            return ResponseEntity.badRequest()
                    .body(createErrorResponse("Invalid shape data: " + e.getMessage()));
        }
    }

    @Operation(summary = "Health check",
              description = "Checks if the service is running properly")
    @GetMapping("/health")
    public ResponseEntity<Map<String, Object>> healthCheck() {
        Map<String, Object> response = new HashMap<>();
        response.put("status", "healthy");
        response.put("service", "Shape Overlap Detection API");
        response.put("timestamp", System.currentTimeMillis());
        return ResponseEntity.ok(response);
    }

    /**
     * Helper method to create error response
     */
    private Map<String, Object> createErrorResponse(String message) {
        Map<String, Object> response = new HashMap<>();
        response.put("success", false);
        response.put("error", message);
        response.put("timestamp", System.currentTimeMillis());
        return response;
    }

    /**
     * Helper method to get shape information
     */
    private Map<String, Object> getShapeInfo(Shape shape) {
        Map<String, Object> info = new HashMap<>();
        info.put("type", shape.getShapeType());
        info.put("x", shape.getX());
        info.put("y", shape.getY());
        info.put("area", shape.getArea());
        return info;
    }
} 