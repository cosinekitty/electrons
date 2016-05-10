/*
    electrons.ts  -  Don Cross  -  http://cosinekitty.com
*/

/// <reference path="jquery.d.ts" />

module Electrons {
    function RadiansFromDegrees(d:number):number {
        return d * (Math.PI / 180);
    }

    class Vector {
        private x:number;
        private y:number;
        private z:number;

        public static Zero:Vector = new Vector(0, 0, 0);

        constructor(_x:number, _y:number, _z:number) {
            this.x = _x;
            this.y = _y;
            this.z = _z;
        }

        public getX():number { return this.x; }
        public getY():number { return this.y; }
        public getZ():number { return this.z; }

        public static Distance(a:Vector, b:Vector):number {
            var dx:number = b.getX() - a.getX();
            var dy:number = b.getY() - a.getY();
            var dz:number = b.getZ() - a.getZ();
            return Math.sqrt(dx*dx + dy*dy + dz*dz);
        }

        public static Cross(a:Vector, b:Vector):Vector {
            return new Vector(
                a.y*b.z - a.z*b.y,
                a.z*b.x - a.x*b.z,
                a.x*b.y - a.y*b.x
            );
        }

        public static Dot(a:Vector, b:Vector):number {
            return a.x*b.x + a.y*b.y + a.z*b.z;
        }

        public dot(other:Vector):number {
            return Vector.Dot(this, other);
        }

        public UnitVector():Vector {
            var mag:number = this.abs();
            return new Vector(this.x/mag, this.y/mag, this.z/mag);
        }

        public absSquared():number {
            return Vector.Dot(this, this);
        }

        public abs():number {
            return Math.sqrt(this.absSquared());
        }

        public add(other:Vector):Vector {
            return new Vector(
                this.x + other.x,
                this.y + other.y,
                this.z + other.z
            );
        }

        public sub(other:Vector):Vector {
            return new Vector(
                this.x - other.x,
                this.y - other.y,
                this.z - other.z
            );
        }

        public neg():Vector {
            return new Vector(-this.x, -this.y, -this.z);
        }

        public mul(scalar:number) {
            return new Vector(scalar*this.x, scalar*this.y, scalar*this.z);
        }

        public RotateX(a:number, b:number):Vector {
            return new Vector(this.x, a*this.y - b*this.z, b*this.y + a*this.z);
        }

        public RotateY(a:number, b:number):Vector {
            return new Vector(a*this.x + b*this.z, this.y, a*this.z - b*this.x);
        }

        public RotateZ(a:number, b:number):Vector {
            return new Vector(a*this.x - b*this.y, b*this.x - a*this.y, this.z);
        }
    }

    class RotationMatrix {
        private constructor(private r:Vector, private s:Vector, private t:Vector) {
        }

        public static Unrotated:RotationMatrix = new RotationMatrix(
            new Vector(1, 0, 0),
            new Vector(0, 1, 0),
            new Vector(0, 0, 1)
        );

        public RotateX(radians:number):RotationMatrix {
            let a:number = Math.cos(radians);
            let b:number = Math.sin(radians);
            return new RotationMatrix(
                this.r.RotateX(a, b),
                this.s.RotateX(a, b),
                this.t.RotateX(a, b));
        }

        public RotateY(radians:number):RotationMatrix {
            let a:number = Math.cos(radians);
            let b:number = Math.sin(radians);
            return new RotationMatrix(
                this.r.RotateY(a, b),
                this.s.RotateY(a, b),
                this.t.RotateY(a, b));
        }

        public RotateZ(radians:number):RotationMatrix {
            let a:number = Math.cos(radians);
            let b:number = Math.sin(radians);
            return new RotationMatrix(
                this.r.RotateZ(a, b),
                this.s.RotateZ(a, b),
                this.t.RotateZ(a, b));
        }

        public Rotate(v:Vector):Vector {
            return new Vector(v.dot(this.r), v.dot(this.s), v.dot(this.t));
        }
    }

    class CameraCoords {
        private hor:number;
        private ver:number;

        constructor (h:number, v:number) {
            this.hor = h;
            this.ver = v;
        }

        public getHor():number { return this.hor; }
        public getVer():number { return this.ver; }
    }

    class Display {
        private rotmat:RotationMatrix = RotationMatrix.Unrotated;

        public constructor(
            private pixelsWide:number,          // number of pixels wide in canvas
            private pixelsHigh:number,          // number of pixels high in canvas
            private zoomFactor:number,          // large values zoom in, small values zoom out
            private parallaxDistance:number)    // for scaling effect of distance on perspective
        {
        }

        public RotateX(radians:number):void {
            this.rotmat = this.rotmat.RotateX(radians);
        }

        public RotateY(radians:number):void {
            this.rotmat = this.rotmat.RotateY(radians);
        }

        public RotateZ(radians:number):void {
            this.rotmat = this.rotmat.RotateZ(radians);
        }

        public Erase(context:CanvasRenderingContext2D):void {
            context.clearRect(0, 0, this.pixelsWide, this.pixelsHigh);
        }

        private Scale(z:number):number {
            return (this.pixelsWide * this.zoomFactor) / (this.parallaxDistance - z);
        }

        private GetCameraCoords(point:Vector):CameraCoords {
            var rot:Vector = this.rotmat.Rotate(point);
            var scale:number = this.Scale(rot.getZ());
            var h:number = scale * (1 + rot.getX());
            var v:number = scale * (1 - rot.getY());
            return new CameraCoords(h, v);
        }

        public DrawSphere(
            context:CanvasRenderingContext2D,
            center:Vector,
            radius:number,
            color:string):void
        {
            let crot:Vector = this.rotmat.Rotate(center);
            let origin:CameraCoords = this.GetCameraCoords(center);
            let theta:number = Math.asin(radius/(this.parallaxDistance - crot.getZ()));
            let cradius:number = this.parallaxDistance * this.Scale(crot.getZ()) * Math.tan(theta);

            context.beginPath();
            context.arc(origin.getHor(), origin.getVer(), cradius, 0, 2*Math.PI, true);
            context.strokeStyle = color;
            context.lineWidth = 1;
            context.stroke();
        }

        public DrawLine(
            context:CanvasRenderingContext2D,
            startpoint:Vector,
            endpoint:Vector,
            color:string):void
        {
            let startcam:CameraCoords = this.GetCameraCoords(startpoint);
            let endcam:CameraCoords = this.GetCameraCoords(endpoint);

            context.beginPath();
            context.moveTo(startcam.getHor(), startcam.getVer());
            context.lineTo(endcam.getHor(), endcam.getVer());
            context.strokeStyle = color;
            context.lineWidth = 1;
            context.stroke();
        }

        public DrawAxis(
            context:CanvasRenderingContext2D,
            name:string,
            startpoint:Vector,
            endpoint:Vector,
            color:string):void
        {
            this.DrawLine(context, startpoint, endpoint, color);
            // FIXFIXFIX - render the name of the axis
        }
    }

    class Particle {
        private force:Vector;
        private position:Vector;

        public constructor(position:Vector) {
            this.position = position.UnitVector();
        }

        public GetPosition():Vector {
            return this.position;
        }

        public ResetForce():void {
            this.force = Vector.Zero;
        }

        public AddForce(other:Particle):void {
            // Force of electrically charged particles:
            // F = k*q1*q2/r^2.
            // We simplify the problem as:  F = 1/r^2.
            // Force is along the direction of the line passing through both.
            let dp:Vector = this.position.sub(other.position);
            let forcemag:number = 1.0 / dp.absSquared();
            let force:Vector = dp.UnitVector().mul(forcemag);
            this.force = this.force.add(force);
            other.force = other.force.sub(force);
        }

        public Migrate(dt:number):void {
            // Update the position of the particle.
            // Move the particle in the direction of the force, but constrain
            // to the surface of the unit sphere.
            // Calculate the tangential component of the force along
            // the sphere's surface.
            // Force = Force_radial + Force_tangential
            // Force_tangential = Force - Force_radial
            // So calculate radial component using dot product and subtract
            // to get tangential component.
            let forceRadial:Vector = this.position.mul(this.force.dot(this.position));
            let forceTangential:Vector = this.force.sub(forceRadial);
            this.position = this.position.add(forceTangential.mul(dt)).UnitVector();
        }
    }

    class Simulation {
        private particleList:Particle[];
        private sphereCenter:Vector = new Vector(0, 0, 0);
        private sphereRadius:number = 1.0;
        private xAxis:Vector = new Vector(1, 0, 0);
        private yAxis:Vector = new Vector(0, 1, 0);
        private zAxis:Vector = new Vector(0, 0, 1);

        public constructor() {
            this.particleList = [];
        }

        public InsertParticle(p:Particle):void {
            this.particleList.push(p);
        }

        public Update(dt:number):void {
            for (let i=0; i < this.particleList.length; ++i) {
                this.particleList[i].ResetForce();
            }

            for (let i=0; i < this.particleList.length - 1; ++i) {
                for (let j=i+1; j < this.particleList.length; ++j) {
                    this.particleList[i].AddForce(this.particleList[j]);
                }
            }

            for (let i=0; i < this.particleList.length; ++i) {
                this.particleList[i].Migrate(dt);
            }
        }

        public Render(display:Display):void {
            var context:CanvasRenderingContext2D = canvas.getContext('2d');
            display.Erase(context);
            display.DrawSphere(context, this.sphereCenter, this.sphereRadius, '#77f');
            display.DrawAxis(context, 'x', this.sphereCenter, this.xAxis, '#f00');
            display.DrawAxis(context, 'y', this.sphereCenter, this.yAxis, '#f00');
            display.DrawAxis(context, 'z', this.sphereCenter, this.zAxis, '#f00');
            for (let i:number = 0; i < this.particleList.length; ++i) {
                display.DrawSphere(context, this.particleList[i].GetPosition(), 0.01, '#000');
            }
        }
    }

    var canvas:HTMLCanvasElement;
    var sim:Simulation;
    var display:Display;
    const FrameDelayMillis:number = 30;
    const SimTimeIncrementSeconds:number = 0.01;
    const ZoomFactor:number = 4.9;
    const ParallaxDistance:number = 10.0;

    function AnimationFrame():void {
        sim.Render(display);
        sim.Update(SimTimeIncrementSeconds);
        display.RotateY(RadiansFromDegrees(0.1));
        window.setTimeout(AnimationFrame, FrameDelayMillis);
    }

    $(document).ready(function(){
        canvas = <HTMLCanvasElement> document.getElementById('SimCanvas');
        sim = new Simulation();
        for (let i:number = 0; i < 17; ++i) {
            let x:number = Math.random();
            let y:number = Math.random();
            let z:number = Math.random();
            sim.InsertParticle(new Particle(new Vector(x, y, z)));
        }
        display = new Display(canvas.width, canvas.height, ZoomFactor, ParallaxDistance);
        display.RotateX(RadiansFromDegrees(-15));
        AnimationFrame();
    });
}
