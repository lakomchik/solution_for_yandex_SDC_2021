  
#include <list>
#include <algorithm>
#include <iostream>
#include <cstddef>
#include <fstream>




class Point {
public:
    Point( int a = 1, int b = 1 ) { x = a; y = b; }
    bool operator ==( const Point& o ) { return o.x == x && o.y == y; }
    Point operator +( const Point& o ) { return Point( o.x + x, o.y + y ); }
    Point operator - (const Point & o) { return Point (x-o.x, y - o.y); }
    void set(const Point & o)
    {
        x = o.x;
        y = o.y;
    }

    void print(){
        std::cout<<" x = "<<x<<" y= "<<y<<std::endl;
    }
    int x, y;
};



bool is_point_exist_in_list(Point point, std::list<Point> points)
{
    bool result = false;
    for (std::list<Point>::iterator i = points.begin(); i!= points.end(); ++i)
    {
        if( point == (*i))
        {
            result = true;
            break;
        }
    }
    return result;
}



void path_print(std::list<Point> path)
{
    for( std::list<Point>::iterator i = path.begin(); i != path.end(); i++ ) {
        std::cout<< "(" << ( *i ).x << ", " << ( *i ).y << ") ";
    }
    std::cout<<std::endl;


}


class Map {
public:
    Map()
    {

    }
    Map(size_t N){

        w = N;
        h = N;
        for (size_t i =0; i < N; ++i)
        {
            for (size_t j = 0; j < N; ++j)
            {
            char c;
                
                std::cin>>c;
                if (c == '#')
                {
                    m[j][i] = '#';
                }
                else
                {
                    m[j][i] = '.';
                }

            }
        
        }
    }

    void print(){

        for (int i = 0; i < h; ++i)
        {
            for (int j =0; j <w; ++j)
            {
                std::cout<<m[j][i];

            }
            std::cout<<std::endl;
        }
    }
   

    int check_empty_cells()
    {
        int free_cells =0;
        for (int i = 0; i < w; ++i)
        {
            for (int j = 0; j< h; ++j)
            {
                if( m[i][j] == '.')
                {
                    ++free_cells;
                }
            }
        }
        return free_cells;
    }
    void map_segmetation(int &segment_size, int nrn, std::list<Point> & starting_points)
    {
        bool is_segments_enough = false;
        char m[w][h];
        while(!is_segments_enough)
        {
            
            for(int i = 0; i <w; ++i)
            {
                for (int j =0; j<w; ++j)
                {
                m[i][j] = this->m[i][j];  
                }
            }

            Point borders[4];
            borders[0] = Point(-1,0);
            borders[1] = Point(0, -1);
            borders[2] = Point(1,0);
            borders[3] = Point (0,1);
            int border_num = 4;
            int group = 1;
            int segment_amount = 0;
            starting_points.clear();
            for (int row = 0; row < w; ++row)
            {
                for (int col = 0; col < w; ++col)
                {
                    int res_x = 0 , res_y = 0;

                    if(m[row][col] == '.')
                    {   
                        Point closed_point(row,col);
                        std::list<Point> opened_points;
                        std::list<Point> closed_points;
                        opened_points.push_back(closed_point);
                        std::list<int> cost;
                        int curr_size = 0;
                        cost.push_back(0);
                        char trace = group%10 + '0';
                      
                        while(opened_points.size() != 0 && curr_size < segment_size)
                        {
                            int min_idx = 0;
                            int min_cost;
                            std::list<int>::iterator x = cost.begin();
                            min_cost = (*x);
                            for (int i = 0; i < cost.size(); ++i)
                            {
                                std::list<int>::iterator l = cost.begin();
                                std::advance(l,i);
                                if( min_cost > (*l))
                                {
                                    min_cost = (*l);
                                    min_idx = i;
                                    
                                }
                                
                            }
                            
                            std::list<Point>::iterator l = opened_points.begin();
                            std::list<int>::iterator y = cost.begin();
                            std::advance(l,min_idx);
                            std::advance(y, min_idx);
                           
                            for (int i = 0; i < border_num; ++i)
                            {
                                Point pretendent  = (*l) + borders[i];
                                if( pretendent.x > -1 && pretendent.y > -1 && pretendent.x < w && pretendent.y < w)
                                {
                                    
                                    if( m[pretendent.x][pretendent.y] == '.' && !is_point_exist_in_list(pretendent, opened_points))
                                    {
                                        opened_points.push_back(pretendent);
                                        int pretendent_cost = (*y)+1;
                                        cost.push_back(pretendent_cost);
                                    }
                                }
                            }
                        
                            m[(*l).x][(*l).y] = trace;
                            closed_points.push_back((*l));
                            opened_points.erase(l);
                            cost.erase(y);
                            ++curr_size;
                        }
                        ++group;
                        std::list<Point>::iterator y = closed_points.begin();
                        std::advance(y, closed_points.size()/2);
                        res_x = (*y).x;
                        res_y = (*y).y;
                        closed_points.clear();
                        char s = 's';
                        m[res_x][res_y] = s;
                        Point starting(res_x, res_y);
                        starting_points.push_back(starting);

                    }
                
                }
            }
            
            if((group-1) > nrn)
            {
                is_segments_enough = true;
            }
            else
            {
                is_segments_enough = true;
            }

            }
      
 


    }

    int operator() ( int x, int y ) { return m[x-1][y-1]; }
    char m[1200][1200];
    int w, h;
};


 
class node {
public:
    bool operator == (const node& o ) { return pos == o.pos; }
    bool operator == (const Point& o ) { return pos == o; }
    bool operator < (const node& o ) { return dist + cost < o.dist + o.cost; }
    Point pos, parent;
    int dist, cost;
};
 
class aStar {
public:
    aStar() {
        neighbours[0] = Point( -1, 0 ); neighbours[1] = Point(  0, 1 );
        neighbours[2] = Point( 1,  0 ); neighbours[3] = Point(  0,  -1 );
    }
 
    int calcDist( Point& p ){
        int x = end.x - p.x, y = end.y - p.y;
        return( x * x + y * y );
    }
 
    bool isValid( Point& p ) {
        return ( p.x >0 && p.y > 0 && p.x < m.w+1 && p.y < m.h+1 );
    }
    void init (Map & m)
    {
        this->m = m;
    }
 
    bool existPoint( Point& p, int cost ) {
        std::list<node>::iterator i;
        i = std::find( closed.begin(), closed.end(), p );
        if( i != closed.end() ) {
            if( ( *i ).cost + ( *i ).dist < cost ) return true;
            else { closed.erase( i ); return false; }
        }
        i = std::find( open.begin(), open.end(), p );
        if( i != open.end() ) {
            if( ( *i ).cost + ( *i ).dist < cost ) return true;
            else { open.erase( i ); return false; }
        }
        return false;
    }
    bool existPoint2(Point & p)
    {
        std::list<node>::iterator i;
        i = std::find( closed.begin(), closed.end(), p );
        if(i != closed.end())
        {
            return true;
        }
        else
        {
            return false;
        }
    }
 
    bool fillOpen( node& n ) {
        int stepCost, nc, dist;
        Point neighbour;
        direction = end - n.pos;
        //std::cout<<"abs x = "<<abs(direction.x)<<"\n";
        //std::cout<<"abs y = "<<abs(direction.y)<<"\n";
        if(abs(direction.x) > abs(direction.y))
        {
            direction.x = direction.x/abs(direction.x);
            direction.y = 0;
        }
        else
        {
            direction.x = 0;
            direction.y = direction.y/abs(direction.y);
        }
        //direction.print();
        neighbour = n.pos + direction;
        if(m(neighbour.x, neighbour.y) !='#' && !existPoint2(neighbour))
        {
            stepCost = 1;
            if( neighbour == end ) return true;
 
            
            nc = stepCost + n.cost;
            dist = calcDist( neighbour );
            if( !existPoint( neighbour, nc + dist ) ) {
                node m;
                m.cost = nc; m.dist = dist;
                m.pos = neighbour;
                m.parent = n.pos;
                open.push_back( m );
               // m.pos.print();
        }
        }
        else
        {
        
            if( neighbour == end ) return true;
            for( int x = 0; x < 4; x++ ) {
                // one can make diagonals have different cost
                //stepCost = x < 4 ? 1 : 1;
                stepCost = 1;
                neighbour = n.pos + neighbours[x];
                if( neighbour == end ) return true;
    
                if( isValid( neighbour ) && m( neighbour.x, neighbour.y ) != '#' ) {
                    nc = stepCost + n.cost;
                    dist = calcDist( neighbour );
                    if( !existPoint( neighbour, nc + dist ) ) {
                        node m;
                        m.cost = nc; m.dist = dist;
                        m.pos = neighbour;
                        m.parent = n.pos;
                        open.push_back( m );
                    }
                }
            }
        }
        return false;
    }
 
    bool search( Point& s, Point& e) {
        if(s == e)
        {
            start =s; end = e;
            //std::cout<<"Point path!\n";
            return true;
        }
        else{
            node n; end = e; start = s;
            n.cost = 0; n.pos = s; n.parent = 0; n.dist = calcDist( s );
            open.push_back( n );
            while( !open.empty() ) {
                node n = open.front();
                open.pop_front();
                closed.push_back( n );
                if( fillOpen( n ) ) return true;
            }
            return false;
        }
    }
 
    int path( std::list<Point>& path ) {
        int cost;
        if(start == end)
        {
            path.push_back(start);
            cost =0;
            //std::cout<<"Point path!\n";
        }
        else
        {
            
            if((abs((start-end).x) +abs((start-end).y))==1)
            {
                path.push_back(start);
                path.push_back(end);
                cost = 1;
            }
            else
            {
                path.push_front( end );
                cost = 1 + closed.back().cost;
                path.push_front( closed.back().pos );
                Point parent = closed.back().parent;
        
                for( std::list<node>::reverse_iterator i = closed.rbegin(); i != closed.rend(); i++ ) {
                    if( ( *i ).pos == parent && !( ( *i ).pos == start ) ) {
                        path.push_front( ( *i ).pos );
                        parent = ( *i ).parent;
                    }
                }
                path.push_front( start );
            }
        }
        return cost;
    }
    void reload(){
        open.clear();
        closed.clear();
    }
 
    Map m; Point end, start, direction;
    Point neighbours[4];
    std::list<node> open;
    std::list<node> closed;
};



Point up(0,-1);
Point left(-1,0);
Point down(0,1);
Point right(1,0);
Point take(0,0);


class Rover{
public:
    Rover(int a, int b){
        position.x = a;
        position.y = b;
        status = false;
        load = false;
    }
    bool is_busy()
    {
        return status;
    }
    Point curr_pos()
    {
        return position;
    }
    std::list<Point> way(){
        return path;
    }
    void append_path(std::list<Point> path){
        for( std::list<Point>::iterator i = path.begin(); i!= path.end(); ++i)
        {
            this->path.push_back(*i);
        }
        status = true;
    }
    char step(){
        char out;
        if(path.size() !=0 )
        {
            if(path.size() !=1 )
            {
                Point current, next;
                std::list<Point>::iterator i = path.begin();
                current = (*i);
                std::advance(i,1);
                next = (*i);
                position = next;
                Point delta = next - current;
                if( delta == up)
                {
                    out = 'U'; 
                }
                if( delta == down)
                {
                   out = 'D'; 
                }
                if( delta == right)
                {
                    out = 'R'; 
                }
                if( delta == left)
                {
                    out = 'L';; 
                }
                if( delta == take && load == false)
                {
                    out = 'T';
                    load = true;
                }
                path.pop_front();
            }
            else
            {
                if(status == true && load == true)
                {
                    out = 'P';
                    path.clear();
                    status = false;
                    load = false;
                }
                else
                {
                    out = 's';
                }
            }
        }
        else
        {
            out = 'S';
        }
        return out;
    }
    std::list<Point> path;
    Point position;
    int idx = -1;
    bool status;
    bool load;
};
 

class Order{
public:
    Order(){

    }
    Order(int xs, int ys, int xe, int ye, int iteration, int & age){
        start.x = xs;
        start.y = ys;
        end.x = xe;
        end.y = ye;
        this->age = age;
        ++ age;
        this->iteration = iteration;

    }
    void set(int xs, int ys, int xe, int ye){
        start.x = xs;
        start.y = ys;
        end.x = xe;
        end.y = ye;

    }
    void print()
    {   
        std::cout<<"Start\n";
        start.print();
        std::cout<<"End\n";
        end.print();
        std::cout<<"Order age: "<<age<<"\n";
    }

    int age;
    int iteration;
    Point start, end;
};











void sort_orders(std::list<Order>& current_orders, std::list<Order>& orders)
{
    std::list<Point> starting_points;
    for(std::list<Order>::iterator i = current_orders.begin(); i != current_orders.end(); ++i)
    {
        starting_points.push_back((*i).start);
    }
    std::list<Order>::iterator i = orders.begin();
    int k=0;
    while(i!=orders.end())
    {
        if (!is_point_exist_in_list((*i).start, starting_points))
        {
            starting_points.push_back((*i).start);
            current_orders.push_back((*i));
            orders.erase(i);
            i = orders.begin();
        }
        else
        {
            ++i;
        }
    }

}


void check_rovers(std::list<Rover>& free_rovers, std::list<Rover>& busy_rovers)
{
    std::list<Rover>::iterator i = busy_rovers.begin();
    while(i!=busy_rovers.end())
    {
        if (!(*i).is_busy())
        {
            free_rovers.push_back(*i);
            busy_rovers.erase(i);
            i = busy_rovers.begin();
        }
        else
        {
            ++i;
        }
    }    
}



void make_path(int& cost, std::list<Point>& path, aStar & as, Point & start, Point & goal)
{
    path.clear();
    as.reload();
    if(as.search(start, goal))
    {
        cost = as.path(path);
    }
    else
    {
        std::cout<<"!!! Error with pathfinding !!!\n";
    }
}



void appoint_rovers(std::list<Rover> & free_rovers,std::list<Rover> & busy_rovers, std::list<Order>& current_orders, aStar & as, bool & is_blocks, int block_start_coord, int block_period)
{
    while(free_rovers.size() != 0 && current_orders.size() != 0)
    {
        int *costes = new int[current_orders.size()];
        std::list<Rover>::iterator j = free_rovers.begin();
        int idx = 0;
         for (std::list<Order>::iterator i = current_orders.begin(); i!= current_orders.end(); ++i)
        {
            int cost;
            std::list<Point> path;
            cost = abs((*j).position.x-(*i).start.x) + abs((*j).position.y-(*i).start.y);
            
            costes[idx] = cost;
            ++ idx;
        }
        
        int min_idx = 0;
        for (int k = 0; k < idx; ++k) 
        {
            if( costes[min_idx] > costes[k])
            {
                min_idx = k;
            }
        }
        
        std::list<Order>::iterator i = current_orders.begin();
        std::advance( i,min_idx);
        std::list<Point> path;
        int cost;
        if(! is_blocks)
        {
            make_path(cost, path, as, (*j).position, (*i).start);
            ((*j).append_path(path));
            make_path(cost, path, as, (*i).start, (*i).end);
            ((*j).append_path(path));
            busy_rovers.push_back((*j));
            free_rovers.erase(j);
            current_orders.erase(i);
        }
        else
        {
            int x_block_start = (((*j).position.x-block_start_coord)/block_period)*block_period + block_start_coord;
           // std::cout<<"The robot\n";
            //(*j).position.print();
           // std::cout<<"x_block_start "<<x_block_start<<"\n"; 
            int y_block_start = (((*j).position.y-block_start_coord)/block_period)*block_period + block_start_coord;
           // std::cout<<"y_block_start "<<y_block_start<<"\n"; 
            Point block_start(x_block_start, y_block_start);
            int x_block_end = (((*i).start.x-block_start_coord)/block_period)*block_period + block_start_coord;
            int y_block_end = (((*i).start.y-block_start_coord)/block_period)*block_period + block_start_coord;
            //std::cout<<"The order start\n";
            //(*i).start.print();
           // std::cout<<"x_block_end "<<y_block_end<<"\n"; 
            //std::cout<<"y_block_end "<<y_block_end<<"\n"; 

            Point block_end(x_block_end, y_block_end);
            Point block_delta(x_block_start, y_block_end);
            //std::cout<<"Starting count this shit\n";
            //block_start.print();
            //block_end.print();

            make_path(cost, path, as, (*j).position, block_start);
            
            //std::cout<<"1st computed\n";
            path.pop_back();
            //path_print(path);
            ((*j).append_path(path));
            
           /* std::cout<<"2nd computed\n";
            path.pop_back();
            path_print(path);
            ((*j).append_path(path));*/
            if(block_start.x != block_end.y && block_start.y != block_end.y)
            {
                make_path(cost, path, as, block_start, block_delta);
                path.pop_back();
                ((*j).append_path(path));
              //  path_print(path);
                make_path(cost, path, as, block_delta, block_end);
                path.pop_back();
                ((*j).append_path(path));
                //path_print(path);
            }
            else
            {
                make_path(cost, path, as, block_start, block_end);
                path.pop_back();
                //path_print(path);
                ((*j).append_path(path));
            }
            
            //std::cout<<"3rd computed\n";
            //block_end.print();
            //(*i).start.print();
            make_path(cost,path,as,block_end, (*i).start);
            (*j).append_path(path);
            //path_print(path);
            /// Then order
            block_start = block_end;
            x_block_end = (((*i).end.x-block_start_coord)/block_period)*block_period + block_start_coord;
            y_block_end = (((*i).end.y-block_start_coord)/block_period)*block_period + block_start_coord;
            block_end.x = x_block_end;
            block_end.y = y_block_end;
            block_delta.x = block_start.x;
            block_delta.y = block_end.y;

            make_path(cost, path, as, (*i).start, block_start);
            path.pop_back();
            //path_print(path);
            ((*j).append_path(path));
            //std::cout<<"4th computed\n";
           // block_start.print();
            //block_end.print();
            if((block_start.x != block_end.y )&&( block_start.y != block_end.y))
            {
               // block_start.print();
               // block_delta.print();
              //  std::cout<<"FUCK\n";
                make_path(cost, path, as, block_start, block_delta);
                path.pop_back();
                ((*j).append_path(path));
                make_path(cost, path, as, block_delta, block_end);
                path.pop_back();
                ((*j).append_path(path));
            }
            else
            {
                //std::cout<<"FUCK2\n";
                make_path(cost, path, as, block_start, block_end);
                path.pop_back();
                //path.pop_back();
                ((*j).append_path(path));
            }
            //std::cout<<"5th computed\n";
            make_path(cost, path, as, block_end, (*i).end);
            (*j).append_path(path);             
            //std::cout<<"6th computed\n";
            //path_print(path);

            
            //path_print((*j).path);



            busy_rovers.push_back((*j));
            free_rovers.erase(j);
            current_orders.erase(i);
        }
    }

}




void robots_step( int robots_amount, std::list<Rover>& free_rovers, std::list<Rover>& busy_rovers, char *msg)
{
    for (int i =0; i<robots_amount;++i)
    {
        msg[i] = 'N';
    }
    for (std::list<Rover>::iterator i = free_rovers.begin(); i !=free_rovers.end(); ++i )
    {
        int k = (*i).idx;
        msg[k] = (*i).step();
    }
    for (std::list<Rover>::iterator i = busy_rovers.begin(); i != busy_rovers.end(); ++i )
    {
        int k = (*i).idx;
        msg[k] = (*i).step();
    }
}

int main( int argc, char* argv[] ) {





    int N, MaxTips, Costc;

    std::cin>>N>>MaxTips>>Costc;
    Map m(N); 
    Map m_copy = m; 
    int k =m.check_empty_cells();

    bool is_blocks = false;
    int block_start = 0;
    int block_period = 0;
    if( N == 180)
    {
        is_blocks = true;
        block_start = 7;
        block_period = 12;
        //std::cout<<"This is blocks\n";

    }
    else
    {
        if( N == 384)
        {
            is_blocks = true;
            block_start = 13;
            block_period = 24;
        }
    
        else
        {
            if (N == 1024)
            {
                is_blocks = true;
                block_start = 17;
                block_period = 32;

            }
        }
    }
    std::list<Point> rovers_starting_points;
	int optimal_robot_area;	
	if (N < 130)
    {
    	optimal_robot_area = 50*50;
    }
	else
    {
    	optimal_robot_area = 30*40;
    }

    int nrn = k/optimal_robot_area;
    if( nrn == 0){
        nrn = 1;
    }
    if ( nrn >100)
    {
        nrn = 100;
    }

    int points_per_area = k/nrn +1;
   
    m_copy.map_segmetation(points_per_area, nrn, rovers_starting_points);
    int iteration_amount, orders_amount, total_orders = 0;
    std::cin>>iteration_amount>>orders_amount;
    aStar as;
	as.init(m);

    std::list<Order> orders;
    std::list<Order> current_orders;
    std::list<Point> orders_start;
    std::list<Rover> free_rovers;
    std::list<Rover> busy_rovers;
    std::list<Point> path;
    int cost;

    if ( nrn > 2 && N>300)    //bozhe praviy, skolko zhe kostiley, stolko nerealizovannih fich, i vse iz zha zhoposhnoi skorosti(((
    {
        nrn = 1;
    }

    Point a1(4,4);
    Point a2(1,4);
    int counter =0;
    std::cout<<nrn<<"\n";
    float returned_size = rovers_starting_points.size();
    float recounter = returned_size/(float)nrn;
    float y = 0;
    
    for (int i = 0; i < nrn ; ++i)
    {
        int idx = y;
        y+=recounter;
        std::list<Point>::iterator l = rovers_starting_points.begin();
        std::advance(l,idx);
        Rover new_rover((*l).x+1,(*l).y+1);
        new_rover.idx = counter;
        free_rovers.push_back(new_rover);
        std::cout<<(*l).y+1<<" "<<(*l).x+1<<"\n";
        ++ counter;
    }
    /*for( float i =0; i < returned_size - recounter; i += recounter)
    {
        //std::cout<<i<<std::endl;
        
        std::list<Point>::iterator l = rovers_starting_points.begin();
        std::advance(l,idx);
        Rover new_rover((*l).x+1,(*l).y+1);
        new_rover.idx = counter;
        free_rovers.push_back(new_rover);
        std::cout<<(*l).y+1<<" "<<(*l).x+1<<"\n";
        ++ counter;
    }*/

   /* for (std::list<Point>::iterator l = rovers_starting_points.begin(); l != rovers_starting_points.end(); ++l)
    {
        Rover new_rover((*l).x+1,(*l).y+1);
        new_rover.idx = counter;
        free_rovers.push_back(new_rover);
        std::cout<<(*l).y+1<<" "<<(*l).x+1<<"\n";
        ++counter;
    }*/








    int rovers_amount = free_rovers.size();


    for (int i =0; i < iteration_amount; ++i)
    {
        int ordersOnIter;
        std::cin>>ordersOnIter;
        for (int j = 0; j < ordersOnIter; ++j)
        {
            int xs, ys, xe, ye;
            std::cin>>ys>>xs>>ye>>xe;
            Order new_order(xs,ys,xe,ye, i, total_orders);
            if( (m(xs,ys) != '#') &&( m(xe,ye) != '#'))
            {
             orders.push_back(new_order);
            }
        }

        sort_orders(current_orders, orders);
        char *msg = new char [rovers_amount*60+1];
        char *time_msg = new char [rovers_amount+1];
        for (int f =0; f < 60; ++f)
        {
            check_rovers(free_rovers,busy_rovers);
            sort_orders(current_orders, orders);
            appoint_rovers(free_rovers, busy_rovers, current_orders, as, is_blocks, block_start, block_period);
            
            
            robots_step(rovers_amount, free_rovers, busy_rovers, time_msg);
            for (int j = 0; j < rovers_amount;++j)
            {
                msg[rovers_amount*f +j] = time_msg[j];
            }
        } 
        for(int k =0 ; k < rovers_amount; ++k)
        {
            for(int u = 0; u <60; ++u)
            {
                std::cout<<msg[k+u*rovers_amount];
            }
            std::cout<<std::endl;

        }   
        




    }



    return 0;
}

